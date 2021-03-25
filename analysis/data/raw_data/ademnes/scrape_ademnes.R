library(tidyverse)
library(rvest)
library(janitor)
library(glue)
library(here)

# Get URLs of site pages from index
read_html("https://www.ademnes.de/db/sites.php") %>%
  html_elements(".textlink_ohne") %>%
  html_attr("href") %>%
  tibble(url = .) %>%
  mutate(url = glue("https://www.ademnes.de/db/{url}")) ->
  ademnes_sites

# For each site we want to retrieve 3–4 tables:
# * Site details
# * Site phases
# * Flora data (if it's there)
# * Fauna data (if it's there)
# As well as a list of references for the flora and/or fauna data.
scrape_ademnes <- function(url) {
  # Get page
  rlang::inform(glue("Scraping {url} ..."))
  page <- read_html(url)

  # Extract elements
  headings <- html_text2(html_elements(page, "h2"))
  tables <- html_elements(page, "table.tabelle")
  reflists <- html_elements(page, "ul")

  # Parse tables
  tables <- html_table(tables, header = TRUE, na.strings = c("", "-"),
                       convert = TRUE)
  names(tables) <- headings

  # Site name
  site_name <- str_remove(headings[1], "Sites => ")

  # Site details
  tables[[1]] %>%
    magrittr::set_colnames(c("X1", "X2")) %>%
    remove_empty("rows") %>%
    pivot_wider(names_from = X1, values_from = X2) %>%
    separate(`ID , Code`, c("id", "site_code"), " , ", convert = TRUE) %>%
    separate(`Latitude, Longitude`, c("latitude", "longitude"), "[^\\d\\.]+",
             convert = TRUE) %>%
    separate(`Date earliest / latest`, c("date_earliest", "date_latest"), " ?/ ?") %>%
    rename(location = "Location",
           geographic_unit = "Geographic Unit",
           recovery_methods = "Recovery Methods",
           preservation_type = "Type of Preservation",
           n_samples = "Nmb. of samples") %>%
    mutate(
      location = na_if(location, "0"),
      latitude = na_if(latitude, 0),
      longitude = na_if(longitude, 0),
      date_earliest = na_if(date_earliest, ""),
      date_latest = na_if(date_latest, "")
    ) ->
    site_details

  # Site phases
  tables[[2]] %>%
    select(phase_code = `Phase Code`,
           phase = `Chronological Specification`,
           period = `Periods`) %>%
    remove_empty("rows") %>%
    extract(period, c("period", "date_earliest", "date_latest"),
            "^(.+?) \\((\\d+?) - (\\d+?)\\)$", convert = TRUE) ->
    phases

  # Flora data
  if ("Flora" %in% headings) {
    flora <- tables[["Flora"]]

    if (ncol(flora) > 6) { # Ah! A nested table.
      flora <- flora[,1:6]
      flora <- filter(flora, `Phase Code` %in% phases$phase_code)
    }

    flora %>%
      select(phase_code = "Phase Code",
             taxon = "Taxon",
             plant_part = "plant part",
             identification_level = "Level of\n      Identific.",
             n = "Nmb.\n      p. Ph.",
             prop = "Prop.\n      p.Ph.") %>%
      remove_empty("rows") %>%
      mutate(across(c(n, prop) & where(is_character),
                    ~parse_number(., locale = locale(decimal_mark = ",")))) ->
      flora
  }
  else flora <- NA

  # Fauna data
  if ("Fauna" %in% headings) {
    fauna <- tables[["Fauna"]]

    if (ncol(fauna) > 9) { # Ah! A nested table.
      fauna <- fauna[,1:9]
      fauna <- filter(fauna, `Phase Code` %in% phases$phase_code)
    }

    fauna %>%
      select(phase_code = "Phase Code",
             taxon = "Taxon",
             n = "Numb.\n      p. Ph.",
             prop = "Prop.\n      p. Ph.",
             weight = "Weight p. Ph.",
             weight_prop = "Weight\n    prop.",
             domestic_prop = "Dom.\n      Prop.",
             domestic_weight_prop = "dom. weight\n    prop.",
             domestic_status = "dom. / wild") %>%
      remove_empty("rows") %>%
      mutate(across(c(n, prop, weight, weight_prop, domestic_prop,
                      domestic_weight_prop) & where(is_character),
                    ~parse_number(., locale = locale(decimal_mark = ",")))) ->
      fauna
  }
  else fauna <- NA

  # References
  references <- html_text2(reflists)
  if (length(references) < 1) references <- NA

  # Return as a single row
  tibble(
    site_name = site_name,
    site_details,
    phases = list(phases),
    flora = list(flora),
    fauna = list(fauna),
    references = list(references)
  )
}

# Scrape all sites (takes a while!)
ademnes <- map_dfr(ademnes_sites$url, scrape_ademnes)

# The same request broken into chunks, for easier debugging
# ademnes1 <- map_dfr(ademnes_sites$url[1:50], scrape_ademnes)
# ademnes2 <- map_dfr(ademnes_sites$url[51:100], scrape_ademnes)
# ademnes3 <- map_dfr(ademnes_sites$url[101:150], scrape_ademnes)
# ademnes4 <- map_dfr(ademnes_sites$url[151:200], scrape_ademnes)
# ademnes5 <- map_dfr(ademnes_sites$url[201:250], scrape_ademnes)
# ademnes6 <- map_dfr(ademnes_sites$url[251:300], scrape_ademnes)
# ademnes7 <- map_dfr(ademnes_sites$url[301:350], scrape_ademnes)
# ademnes8 <- map_dfr(ademnes_sites$url[351:400], scrape_ademnes)
# ademnes9 <- map_dfr(ademnes_sites$url[401:length(ademnes$url)], scrape_ademnes)
# ademnes <- bind_rows(ademnes1, ademnes2, ademnes3, ademnes4, ademnes5, ademnes6,
#                      ademnes7, ademnes8, ademnes9)

# Serialise .RData and write to TSV
save(ademnes, file = here("analysis/data/raw_data/ademnes/ademnes.RData"))

ademnes_site <- select(ademnes, id, site_name, site_code:n_samples)
write_tsv(ademnes_site, here("analysis/data/raw_data/ademnes/ademnes_site.tsv"))

select(ademnes, id, site_name, site_code, phases) %>%
  unnest(phases) ->
  ademnes_phase
write_tsv(ademnes_phase, here("analysis/data/raw_data/ademnes/ademnes_phase.tsv"))

select(ademnes, id, site_name, site_code, flora) %>%
  mutate(
    flora = map(flora, ~ if(!all(is.na(.x))) {
      mutate(.x, identification_level = as.character(identification_level))
    }
    else NA
    )
  ) %>%
  unnest(flora) %>%
  select(-flora) %>%
  drop_na(phase_code) ->
  ademnes_flora
write_tsv(ademnes_flora, here("analysis/data/raw_data/ademnes/ademnes_flora.tsv"))

select(ademnes, id, site_name, site_code, fauna) %>%
  unnest(fauna) %>%
  select(-fauna) %>%
  drop_na(phase_code) ->
  ademnes_fauna
write_tsv(ademnes_fauna, here("analysis/data/raw_data/ademnes/ademnes_fauna.tsv"))

select(ademnes, id, site_name, site_code, references) %>%
  unnest(references) %>%
  rename(reference = "references") %>%
  drop_na(reference) ->
  ademnes_reference
write_tsv(ademnes_reference, here("analysis/data/raw_data/ademnes/ademnes_reference.tsv"))
