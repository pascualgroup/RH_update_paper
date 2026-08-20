  surat_temp <- analyze_profile(
    # "aggregated_results/profile_Surat_Max_Temp_1997_2016_5_9.csv",
    "aggregated_results/all_profiles.csv",
    city = "Surat", covariate_label = "Maximum temperature", output_prefix = "Surat_Max_Temp"
  )

    read_csv( "aggregated_results/profile_Surat_Max_Temp_1997_2016_5_9.csv") |>
    filter(bH <= 0) |>
    write_csv(aggregated_results/profile_Surat_Max_Temp_1997_2016_5_9.csv)

  print(plot_profile_group(surat_temp, base_size = 12, label_size = 5))

  print(plot_profile_group_transformed(surat_temp, base_size = 12, label_size = 5))

  save_profile_plot(surat_temp, param_groups$transmission, "figures/profile_Surat_Max_Temp_transmission.pdf")
  save_profile_plot(surat_temp, param_groups$initial, "figures/profile_Surat_Max_Temp_initial.pdf")
  save_profile_plot(surat_temp, param_groups$immunity, "figures/profile_Surat_Max_Temp_immunity.pdf")


