# Draft table for SI2
flora |>
  drop_na(founder_crop) |>
  mutate(plant_part = if_else(str_detect(taxon_source, "grain"), "grain", "undifferentiated")) |>
  group_by(founder_crop, period, site_name, phase_code, plant_part) |>
  summarise(
    age = first(paste0(age_start, "–", age_end, " BP")),
    total_remains = sum(n),
    references = paste(reference, collapse = ";")
  ) |>
  mutate(
    references = map(str_split(references, ";"), unique),
    references = map_chr(references, paste, collapse = "; ")
  )|>
  write_csv("analysis/figures/SI2_table.csv")

