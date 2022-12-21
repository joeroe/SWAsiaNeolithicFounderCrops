# Draft table for SI2
flora |>
  drop_na(founder_crop) |>
  group_by(founder_crop, period, site_name, phase_code) |>
  summarise(
    age = first(paste0(age_start, "–", age_end, " BP")),
    total_remains = sum(n),
  ) |>
  write_csv("analysis/figures/SI2_table.csv")

