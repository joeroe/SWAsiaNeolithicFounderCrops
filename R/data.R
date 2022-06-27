#' @export
collate_flora <- function(origins_path, ademnes_path, compag_path, taxa_path,
                          region, start_bp, end_bp, quiet = TRUE) {
  # READ DATA
  origins <- read_origins(origins_path)
  ademnes <- read_ademnes(ademnes_path)
  compag <- read_compag2018(compag_path)

  # EXCLUDE NON-SITES
  ademnes <- filter(ademnes, !site_name %in% c("unknown1", "Asvan region"))

  # FILTER BY REGION
  # NB. Not necessary for COMPAG
  origins <- sf::st_as_sf(origins, coords = c("longitude", "latitude"),
                          remove = FALSE)
  sf::st_crs(origins) <- 4326
  origins <- origins[region,]
  origins <- tibble::as_tibble(origins)

  ademnes <- filter(ademnes, location != "Greece" | is.na(location))
  ademnes <- sf::st_as_sf(ademnes, coords = c("longitude", "latitude"),
                          remove = FALSE)
  sf::st_crs(ademnes) <- 4326
  ademnes <- ademnes[region,]
  ademnes <- tibble::as_tibble(ademnes)

  # NORMALISE
  # Dates to ages BP
  origins <- mutate(origins,
                    age_start = -start_bc + 1950,
                    age_end = -end_bc + 1950)
  ademnes <- mutate(ademnes,
                    age_start = date_earliest + 1950,
                    age_end = date_latest + 1950)
  compag <- mutate(compag,
                   age_start = -date_start + 1950,
                   age_end = -date_end + 1950)

  # Aggregate ORIGINS samples by phase
  origins %>%
    group_by(site_code, phase_code, taxon) %>%
    summarise(
      id = glue("{min(id)}–{max(id)}"),
      site = first(site), site_code = first(site_code),
      latitude = first(latitude), longitude = first(longitude),
      phase = first(phase), phase_description = first(phase_description),
      samples = glue_collapse(sample, sep = ";"),
      n = sum(n, na.rm = TRUE), mnpp = sum(mnpp, na.rm = TRUE),
      age_start = max(age_start), age_end = min(age_end),
      family = first(family), genus = first(genus),
      reference = glue_collapse(reference, ";"),
      .groups = "drop"
    ) ->
    origins

  # N to proportions
  origins %>%
    dplyr::group_by(phase_code) %>%
    dplyr::mutate(prop = n / sum(n, na.rm = TRUE)) %>%
    ungroup() ->
    origins
  compag %>%
    dplyr::group_by(phase_code) %>%
    dplyr::mutate(prop = n / sum(n, na.rm = TRUE)) %>%
    ungroup() ->
    compag

  # Standardise columns
  origins <- dplyr::transmute(
    origins,
    source = "ORIGINS",
    source_id = as.character(id),
    source_site_name = site,
    site_name = swap_control_site_name(source_site_name, quiet = TRUE),
    latitude,
    longitude,
    phase,
    phase_description,
    phase_code,
    age_start,
    age_end,
    family,
    genus,
    taxon,
    n,
    prop,
    reference
  )

  ademnes <- dplyr::transmute(
    ademnes,
    source = "ADEMNES",
    source_id = as.character(id),
    source_site_name = site_name,
    site_name = swap_control_site_name(source_site_name, quiet = TRUE),
    latitude,
    longitude,
    phase,
    phase_description = period,
    phase_code,
    age_start,
    age_end,
    family = NA,
    genus = NA,
    taxon,
    n,
    prop,
    reference = NA
  )

  compag <- dplyr::transmute(
    compag,
    source = "COMPAG",
    source_id = NA_character_,
    source_site_name = site_abot,
    site_name = site_canon,
    latitude,
    longitude,
    phase = NA,
    phase_description = NA,
    phase_code,
    age_start,
    age_end,
    family = NA,
    genus = NA,
    taxon,
    n,
    prop,
    reference = references
  )

  # FILTER DUPLICATES
  # Preferring the database with the greatest detail,
  # i.e. ORIGINS > COMPAG > ADEMNES
  compag <- filter(compag, !site_name %in% origins$site_name)
  ademnes <- filter(ademnes,
                    !site_name %in% origins$site_name,
                    !site_name %in% compag$site_name)

  # BIND
  flora <- dplyr::bind_rows(origins, ademnes, compag)

  # FILTER BY PERIOD
  flora <- dplyr::filter(flora, age_end <= start_bp & age_start >= end_bp)


  # STANDARDISE AND CLASSIFY TAXA
  taxa <- read_tsv(taxa_path, col_types = "ccccccccccl")
  flora <- select(flora, -family, -genus)
  flora <- left_join(flora, taxa, by = c("taxon" = "variant"))
  flora <- rename(flora, taxon_source = taxon, taxon_detail = canon, taxon = group)

  return(flora)
}

#' @export
read_origins <- function(path) {
  # Read tables
  sites <- readr::read_csv(fs::path(path, "1_Site.txt"),
                           col_types = "cccdd",
                           locale = readr::locale(encoding = "latin1"))
  phases <- readr::read_csv(fs::path(path, "2_SitePhase.txt"),
                            col_types = "ccccddc",
                            locale = readr::locale(encoding = "latin1"))
  samples <- readr::read_csv(fs::path(path, "3_Records-Samples.txt"),
                             col_types = "ccccccccccdccc")
  references <- readr::read_csv(fs::path(path, "4_References.txt"),
                                col_types = "ccccccccicc")
  taxa <- readr::read_csv(fs::path(path, "5_Taxa_Entries.txt"),
                          col_types = "ciccccclcdci")
  plant_parts <- readr::read_csv(fs::path(path, "6_PlantPart_Categories.txt"),
                                 col_types = "cc")
  taxa_standard <- readr::read_csv(fs::path(path, "7_Taxa_Standardisation.txt"),
                                   col_types = "cc")
  taxa_groups <- readr::read_csv(fs::path(path, "8_Taxa_Groups.txt"),
                                 col_types = "ccccc")

  # Join tables
  taxa %>%
    dplyr::left_join(samples, by = "RecordID") %>%
    dplyr::left_join(phases, by = "SitePhase-ID") %>%
    dplyr::left_join(sites, by = c("SiteName (SPh)" = "SiteName")) %>%
    dplyr::left_join(plant_parts, by = "PlantPart") %>%
    dplyr::left_join(taxa_standard, by = "OriginalTaxon") %>%
    dplyr::left_join(taxa_groups, by = "ProjectTaxon") ->
    origins

  # Reorder and rename columns
  origins <- dplyr::select(origins,
                           id = Spcode,
                           site = `SiteName (SPh)`,
                           site_code = SiteCode,
                           region = `Region-Coarse`,
                           latitude = Latitude,
                           longitude = Longitude,
                           phase_code = `SitePhase-ID`,
                           phase = PhaseName,
                           phase_description = PhaseDescription,
                           period = MainPeriod,
                           start_bc = ApproxStartBC,
                           end_bc = ApproxEndBC,
                           sample = RecordID,
                           context = Context,
                           context_detail = ContextDetail,
                           record_type = RecordType,
                           flot = Flot,
                           area = Area,
                           unit = Unit,
                           level = Level,
                           volume = Volume,
                           recovery_method = RecoveryMethod,
                           sample_notes = Notes,
                           exclusion_notes = `Exclusion notes (if any)`,
                           reference = RefCode,
                           family = Family,
                           genus = Genus,
                           taxon = ProjectTaxon,
                           taxon_type = IDtype,
                           domesticated = DomProgWild,
                           original_taxon = OriginalTaxon,
                           identification_level = IDlevel,
                           identification_notes = IDnotes,
                           plant_part = PlantPart_Cat,
                           original_plant_part = PlantPart,
                           preservation = Preservation,
                           quantified = Quantification,
                           quantification_method = ScoringSystem,
                           quantification_notes = QuantificationNotes,
                           n = NumberInd,
                           mnpp = MNPP)

  # Expand coded and unclear values
  origins %>%
    dplyr::mutate(preservation = dplyr::recode(preservation,
                                               c = "charred",
                                               d = "desiccated",
                                               i = "impression",
                                               w = "waterlogged",
                                               m = "mineralised")) %>%
    dplyr::mutate(taxon_type = dplyr::recode(taxon_type,
                                             `Single ID` = "single",
                                             `Multi ID` = "multiple")) ->
    origins

  # Make missing values explicit
  origins %>%
    dplyr::mutate(genus = dplyr::na_if(genus, "Above Genus ID"),
                  taxon = dplyr::na_if(taxon, "Unidentified"),
                  plant_part = dplyr::na_if(plant_part, "Other/Unknown"),
                  domesticated = dplyr::na_if(domesticated, ".Indet")) ->
    origins

  # Interpolate missing data for "CHGO-iX"
  origins$site[origins$phase_code == "CHGO-iX"] <- "Chogha Golan"
  origins$site_code[origins$phase_code == "CHGO-iX"] <- "CHGO"
  origins$region[origins$phase_code == "CHGO-iX"] <- "E Fertile Crescent"
  origins$latitude[origins$phase_code == "CHGO-iX"] <- 33.41
  origins$longitude[origins$phase_code == "CHGO-iX"] <- 46.26
  origins$start_bc[origins$phase_code == "CHGO-iX"] <- -8900
  origins$end_bc[origins$phase_code == "CHGO-iX"] <- -8700

  # TODO: references? to BibTex?

  return(origins)
}

#' @export
read_ademnes <- function(path) {
  site <- readr::read_tsv(fs::path(path, "ademnes_site.tsv"),
                          col_types = "iccccddcccci")
  phase <- readr::read_tsv(fs::path(path, "ademnes_phase.tsv"),
                           col_types = "icccccdd")
  flora <- readr::read_tsv(fs::path(path, "ademnes_flora.tsv"),
                           col_types = "iccccccdd")

  # Normalise column names
  site <- dplyr::rename(site,
                        site_date_earliest = date_earliest,
                        site_date_latest = date_latest)

  # Normalise percentages
  flora <- dplyr::mutate(flora, prop = prop / 100)

  # Join tables
  ademnes <- dplyr::right_join(site, phase,
                               by = c("id", "site_name", "site_code"))
  ademnes <- dplyr::right_join(ademnes, flora,
                               by = c("id", "site_name", "site_code", "phase_code"))

  # Add missing coordinates (from Wikidata)
  ademnes <- dplyr::mutate(
    ademnes,
    latitude = dplyr::recode(site_name,
                             "Dhiban" = 31.50,
                             "Hirbet Iskander" = 31.56,
                             "Jaffa" = 32.05,
                             "Pella" = 32.45,
                             "Sheikh-e Abad" = 34.61,
                             "Sidon" = 33.56,
                             "Tel Beth Yerah" = 32.72,
                             "Tel Malhat" = 31.22,
                             "Tell el-Fukhar" = 32.56,
                             "Tell es-Safi/Gath" = 31.70,
                             "Tell Hadar" = 32.85,
                             "Tell Miqne" = 31.78,
                             "Zahrat adh-Dhra 1" = 31.255,
                             .default = latitude
    ),
    longitude = dplyr::recode(site_name,
                              "Dhiban" = 35.78,
                              "Hirbet Iskander" = 35.77,
                              "Jaffa" = 34.75,
                              "Pella" = 35.61,
                              "Sheikh-e Abad" = 47.27,
                              "Sidon" = 35.40,
                              "Tel Beth Yerah" = 35.57,
                              "Tel Malhat" = 35.03,
                              "Tell el-Fukhar" = 35.85,
                              "Tell es-Safi/Gath" = 34.84,
                              "Tell Hadar" = 35.65,
                              "Tell Miqne" = 34.85,
                              "Zahrat adh-Dhra 1" = 35.57,
                              .default = longitude
    ),
  )


  # Return, dropping empty column n_samples
  dplyr::select(ademnes, -n_samples)
}

#' @export
read_compag2018 <- function(path) {
  # Read tables from CSV
  # Character encoding is unknown and appears to differ between files, so the
  # locales below are a guess based on readr::guess_encoding()
  readr::read_csv(fs::path(path, "AbotData.csv"),
                  col_types = readr::cols(taxon = readr::col_character(),
                                          .default = readr::col_integer()),
                  locale = readr::locale(encoding = "latin1"),
                  n_max = 42) ->
    compag_data
  readr::read_csv(fs::path(path, "AbotSitesrefs.csv"),
                  col_names = c("site", "country", "region", "Latitude",
                                "Longitude", "Georef rank", "Abot sampling",
                                "dating", "Phase", "start date", "end date",
                                "Median age", paste0("X", 1:12), "references"),
                  col_types = "cccddicccdddccccccccccccc",
                  na = c("", "NA", "-"),
                  skip = 1,
                  locale = readr::locale(encoding = "windows-1252")) ->
    compag_site

  # Clean and tidy data table
  compag_data %>%
    dplyr::slice(2:nrow(.)) %>%
    dplyr::select(-TOTALS) %>%
    tidyr::pivot_longer(-tidyselect::any_of("taxon")) %>%
    dplyr::transmute(
      site = stringr::str_extract(name, "[^\\(]+"),
      site = stringr::str_remove(site, " $"),
      phase_code = stringr::str_extract(name, "\\(.+\\)"),
      phase_code = stringr::str_remove_all(phase_code, "[()]"),
      taxon = taxon,
      n = value
    ) ->
    compag_data

  # Clean and tidy site table
  compag_site %>%
    tidyr::unite("characters", X1:X11, sep = ";", na.rm = TRUE) %>%
    dplyr::transmute(
      site,
      country,
      region,
      latitude = Latitude,
      longitude = Longitude,
      georef_rank = `Georef rank`,
      sampling = `Abot sampling`,
      dating,
      phase = Phase,
      date_start = `start date`,
      date_end = `end date`,
      date_mid = `Median age`,
      characters,
      references) ->
    compag_site

  # There are many discrepancies between site names in data table and in the
  # site table, so standardise before joining.
  # And not all sites in the data are present in the site table!
  compag_data$site[compag_data$phase_code == "5960-TSAI"] <- "Sabi Abyad I"
  compag_data$site[compag_data$phase_code == "6300-TSAII"] <- "Sabi Abyad II"
  compag_site$site[compag_site$site == "Sabi Abyad"] <- "Sabi Abyad I"
  compag_site$site[compag_site$site == "Yarim Tepe"] <- "Yarim Tepe I"
  compag_site$site[compag_site$site == "Tell ‘Abr 3"] <- "Tell 'Abr"
  compag_site <- compag_site[compag_site$site != "Tall-e Jari",] # No data?
  compag_site <- dplyr::distinct(compag_site)

  compag_data <- dplyr::mutate(
    compag_data,
    site_canon = swap_control_site_name(site, quiet = TRUE)
  )
  compag_site <- dplyr::mutate(
    compag_site,
    site_canon = swap_control_site_name(site, quiet = TRUE)
  )

  # The data is actually recorded by phase, but there is no explicit link
  # between the phase codes in the data table and the site table. However, we
  # can partly reconstruct the join using a combination of the site name and
  # median age.
  compag_data <- dplyr::mutate(
    compag_data,
    date_mid = stringr::str_extract(phase_code, "[^- ]+"),
    date_mid = suppressWarnings(-as.numeric(date_mid))
  )
  compag <- right_join(compag_site, compag_data, by = c("site_canon", "date_mid"),
                       suffix = c("_site", "_abot"))

  # Clean up final table and return
  dplyr::select(compag,
                site_canon,
                site_abot,
                site_site,
                country,
                region,
                latitude,
                longitude,
                georef_rank,
                sampling,
                phase,
                phase_code,
                dating,
                date_start,
                date_end,
                date_mid,
                taxon,
                n,
                references)
}

#' @export
read_colledge2004 <- function(path, europe_only = FALSE) {
  # Read
  phases <- readr::read_csv(fs::path(path, "phases.csv"), col_types = "cccc")
  samples <- readr::read_csv(fs::path(path, "samples.csv"), col_types = "ccc")
  sites <- readr::read_csv(fs::path(path, "sites.csv"), col_types = "ccc")
  taxa <- readr::read_csv(fs::path(path, "taxa.csv"), col_types = "ccccc")

  # Hot fixes
  samples$PHASE <- recode(
    samples$PHASE,
    "KnH" = "KNH",
    "NHe" = "Nhe",
    "PCa" = "Pca",
    "SGy" = "Sgy",
    "SMa" = "Sma"
  )

  # Join
  colledge2004 <- left_join(samples, taxa, by = "TAXON")
  colledge2004 <- right_join(phases, colledge2004, by = "PHASE")
  colledge2004 <- right_join(sites, colledge2004, by = "SITE")

  # Normalise names
  names(colledge2004) <- tolower(names(colledge2004))

  # Filter if desired
  if (isTRUE(europe_only)) {
    europe <- setdiff(
      unique(colledge2004$country),
      c("Iran", "Iraq", "Israel", "Jordan", "Palestinian Territories", "Syria",
        "Turkey")
    )
    colledge2004 <- dplyr::filter(colledge2004, country %in% europe)
  }

  # Return
  colledge2004
}
