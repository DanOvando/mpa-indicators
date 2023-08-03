



# setup -------------------------------------------------------------------
library(tidyverse)

library(marlin)

library(sf)

library(rerddap)

library(rnaturalearth)

library(here)

library(janitor)

library(progress)

library(furrr)

options(dplyr.summarise.inform = FALSE)

theme_set(theme_minimal(base_size = 12))

run_experiments <- FALSE

results_name <- "v1.0"

results_dir <- file.path("results", results_name)

plot_dir <- file.path(results_dir,"plots")

plot_width <- 8

plot_height = 5

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  
  dir.create(plot_dir, recursive = TRUE)
  
}

set.seed(24)

seasons <- 2

experiment_workers <- 6

dx <-  100000

dy <- 20000

patch_area <- 1

n_states <- 100

resolution <- c(20, 20)

patches <- prod(resolution)

foos <- list.files(here("R"))

walk(foos, ~ source(here("R", .x)))

# load data ---------------------------------------------------------------

# get maps
campas <-
  sf::st_read(here("data", "CA_MPA_boundaries", "ds582"),) |>
  janitor::clean_names()

campas_crs <- sf::st_crs(campas)

ca <-
  rnaturalearth::ne_states(country = "united states of america", returnclass = "sf") |>
  filter(name == "California") |>
  janitor::clean_names() |>
  sf::st_transform(campas_crs)

ca_crs <- sf::st_crs(ca)


ca_raster <- stars::st_rasterize(ca["scalerank"])$scalerank


ca_raster |>
  as.data.frame() |>
  mutate(x = 1:length(V1)) |>
  pivot_longer(
    starts_with("V"),
    names_to = "y",
    values_to = "sst",
    names_prefix = list(y = "V"),
    names_transform = list(y = as.integer)
  ) |>
  ggplot(aes(x,-y, color = sst)) +
  geom_point()


campas |>
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(aes(fill = type))

campas_raster <-
  stars::st_rasterize(campas["acres"], dx = dx, dy = 200)

# get environmental data

# get bathymetry

if (!file.exists(here('data', "bathy.rds"))) {
  a = ed_search(query = 'bathymetry', which = "grid")
  
  bath_dat <- info('srtm30plus_v11_bathy')
  
  (bathy <- griddap(
    bath_dat,
    latitude = c(32, 43),
    longitude = c(-125,-116),
    stride = 20,
    fmt = "csv"
  ))
  
  
  bathy <- bathy %>%
    filter(!is.na(elev) & elev > -2000)
  
  bathy <- st_as_sf(
    bathy,
    coords = c("longitude", "latitude"),
    remove = FALSE,
    crs = "WGS84"
  ) |>
    st_transform(campas_crs)
  
  
  write_rds(bathy, file = here("data", "bathy.rds"))
} else {
  bathy <-  read_rds(file = here("data", "bathy.rds"))
  
  
}


bathy_raster <-
  stars::st_rasterize(bathy["elev"], dx = 100000, dy = 20000)$elev


bathy_coords <- sf::st_coordinates(bathy)

bathy <- bathy |>
  bind_cols(bathy_coords)


bathy |>
  ggplot() +
  geom_sf(aes(color = elev)) +
  geom_sf(data = ca) +
  geom_sf(data = campas) +
  scale_color_viridis_c(option = "magma")


# get SST

if (!file.exists(here('data', "sst.rds"))) {
  sst_dat <- info('jplMURSST41')
  
  sst <- griddap(
    sst_dat,
    latitude = c(32, 43),
    longitude = c(-125, -116),
    stride = 10,
    time = c("2010-01-01", "2011-01-01"),
    fmt = "csv"
  )
  
  write_rds(sst, file = here("data", "sst.rds"))
  
  
} else {
  sst <- read_rds(file = here("data", "sst.rds"))
  
}

sst <- st_as_sf(
  sst,
  coords = c("longitude", "latitude"),
  remove = FALSE,
  crs = "WGS84"
) |>
  st_transform(campas_crs) |>
  mutate(year = lubridate::year(time),
         month = lubridate::month(time))

sst_coords <- sf::st_coordinates(sst)

sst <- sst |>
  bind_cols(sst_coords)

sst_raster <-
  (stars::st_rasterize(sst["analysed_sst"], dx = 100000, dy = 20000)$analysed_sst)

sst_raster |>
  as_tibble() |>
  mutate(x = 1:length(V1)) |>
  pivot_longer(
    starts_with("V"),
    names_to = "y",
    values_to = "sst",
    names_prefix = list(y = "V"),
    names_transform = list(y = as.integer)
  ) |>
  ggplot(aes(x,-y, color = sst)) +
  geom_point() +
  scale_color_viridis_c(option = "magma")

# sst |>
#   filter(!is.na(analysed_sst), year == max(year)) |>
#   group_by(year, month, X, Y) |>
#   summarise(mean_sst = mean(analysed_sst)) |>
#   ggplot() +
#   geom_sf(aes(color = mean_sst)) +
#   geom_sf(data = ca) +
#   geom_sf(data = campas) +
#   scale_color_viridis_c(option = "magma") +
#   facet_wrap( ~ month)


# set up simulations ------------------------------------------------------
if (run_experiments) {
  
campas_bbox <- sf::st_bbox(campas)

lon_range <- seq(campas_bbox$xmin, campas_bbox$xmax, by = 10000)

lat_range <- seq(campas_bbox$ymin, campas_bbox$ymax, by = 10000)

state_experiments <-
  tibble(
    kiss = sample(c(TRUE, FALSE), n_states, replace = TRUE),
    mpa_response = sample(c("stay", "leave"), n_states, replace = TRUE),
    sigma_centroid = runif(n_states, .1 * resolution[1], resolution[1] ^ 2 / 20),
    sigma_hab = runif(n_states, .1 * resolution[1], resolution[1] ^ 2 / 20),
    spatial_q = sample(
      c(FALSE, TRUE),
      n_states,
      replace = TRUE,
      prob = c(3, 1)
    ),
    spatial_allocation = sample(c("ppue", "rpue"), n_states, replace = TRUE),
    fleet_model = sample(c("open access", "constant effort"), n_states, replace = TRUE)
  ) %>%
  mutate(state_id = 1:nrow(.))

critters <-
  tibble(
    scientific_name = c(
      "paralabrax clathratus",
      "ophiodon elongatus",
      "scorpaenichthys marmoratus"
    )
  ) %>%
  mutate(centroid = NA)

create_critter_habitats <-
  function(sigma_centroid,
           sigma_hab = 0.2,
           base_centroid = c(resolution[1] / 2, resolution[1] / 2) ,
           critters,
           resolution) {
    base_layer <-
      tidyr::expand_grid(x = 1:resolution[1], y = 1:resolution[2])
    
    
    centroid_index <-
      which(base_layer$x == base_centroid[1] &
              base_layer$y == base_centroid[2])
    
    
    critters$centroid <-
      pmin(nrow(base_layer), pmax(1, round(
        rnorm(nrow(critters), centroid_index, sigma_centroid)
      )))
    
    critters$habitat <-
      vector(mode = "list", length = nrow(critters))
    
    for (i in 1:nrow(critters)) {
      tmp <- base_layer
      
      tmp <- base_layer %>%
        mutate(c_x = base_layer$x[critters$centroid[i]],
               c_y = base_layer$y[critters$centroid[i]]) %>%
        mutate(distance = sqrt((c_x - x) ^ 2 + (c_y - y) ^ 2)) %>%
        mutate(habitat = resolution[1] * dnorm(distance, 0, sigma_hab)) %>%
        select(x, y, habitat)
      
      tmp$habitat <- tmp$habitat - min(tmp$habitat)
      
      
      critters$habitat[[i]] <-  tmp
      
      
    }
    
    return(critters)
  }

port_locations <-
  tibble(x = c(1, resolution[1]), y = c(1, resolution[2])) # coordinates to test impact of highly disparate ports


# convert maps into grid form for simulation

# some ideas here. get the SST data and then add a bit of buffer to it
# Then, clip the CA layer to that just to have the edges as well

# prepare critters

state_experiments <- state_experiments %>%
  mutate(
    habitats = map2(
      sigma_centroid,
      sigma_hab,
      create_critter_habitats,
      critters = critters,
      resolution = resolution
    )
  ) %>%
  unnest(cols = habitats) %>%
  ungroup() |>
  mutate(
    seasonal_movement = sample(c(FALSE, TRUE), length(state_id), replace = TRUE),
    spawning_aggregation = sample(c(TRUE, FALSE), length(state_id), replace = TRUE),
    spawning_season = sample(1:seasons, length(state_id), replace = TRUE),
    f_v_m = runif(length(state_id), 0.01, 0.2),
    adult_diffusion = runif(length(state_id),
                            min = .25 * resolution,
                            max = patches),
    steepness = runif(length(state_id), min = 0.3, max = 1),
    ssb0 = runif(length(state_id), min = 10, max = 10),
    recruit_diffusion = runif(length(state_id),
                              min = 0,
                              max = patches),
    hyper = runif(length(state_id), 1, 2),
    density_dependence = sample(
      c(
        "global_habitat",
        "local_habitat",
        "pre_dispersal",
        "post_dispersal"
      ),
      length(state_id),
      replace = TRUE
    )
  ) %>%
  mutate(
    spawning_season = ifelse(spawning_aggregation, spawning_season, NA),
    ontogenetic_shift = sample(c(FALSE, TRUE), length(state_id), replace = TRUE)
  )


state_experiments <- state_experiments %>%
  mutate(critter = pmap(
    list(
      sciname = scientific_name,
      habitat = habitat,
      seasonal_movement = seasonal_movement,
      spawning_aggregation = spawning_aggregation,
      spawning_season = spawning_season,
      f_v_m = f_v_m,
      adult_diffusion = adult_diffusion,
      recruit_diffusion = recruit_diffusion,
      density_dependence = density_dependence,
      hyper = hyper,
      ontogenetic_shift = ontogenetic_shift,
      steepness = steepness,
      ssb0 = ssb0,
      kiss = kiss
    ),
    create_experiment_critters,
    resolution = resolution,
    seasons = seasons,
    .progress = TRUE
  ))

# aggregate into lists of fauna

state_experiments <- state_experiments %>%
  group_by(state_id) %>%
  nest() %>%
  mutate(fauna = map(data, ~ .x$critter %>% set_names(.x$scientific_name)))

state_experiments <- state_experiments %>%
  ungroup() %>%
  mutate(use_ports = sample(c(FALSE, TRUE), nrow(.), replace = TRUE)) %>%
  mutate(
    fleet = pmap(
      list(
        fauna = fauna,
        state = data,
        use_ports = use_ports
      ),
      create_experiment_fleet,
      port_locations = port_locations,
      resolution = resolution
    )
  )

# add in starting conditions
init_condit <- function(fauna, fleets, years = 100) {
  starting_trajectory <-
    simmar(fauna = fauna,
           fleets = fleets,
           years = years)
  
  # plot_marlin(check)
  
  starting_conditions <-
    starting_trajectory[(length(starting_trajectory) - seasons + 1):length(starting_trajectory)]
  
  proc_starting_conditions <-
    process_marlin(starting_conditions, keep_age = FALSE)
  
  out <- list(starting_conditions = starting_conditions,
              proc_starting_conditions = proc_starting_conditions)
  
}

state_experiments <- state_experiments %>%
  mutate(tmp = map2(fauna, fleet, init_condit, .progress = TRUE))

state_experiments$starting_conditions <-
  map(state_experiments$tmp, "starting_conditions")

state_experiments$proc_starting_conditions <-
  map(state_experiments$tmp, "proc_starting_conditions")

state_experiments <- state_experiments %>%
  select(-tmp)

placement_experiments <- expand_grid(
  placement_strategy = c("target_fishing", "area"),
  prop_mpa = seq(0, 1, by = 0.05),
  critters_considered = seq(length(state_experiments$fauna[[1]]), length(state_experiments$fauna[[1]]), by = 1),
  placement_error = c(0)
) %>%
  group_by_at(colnames(.)[!colnames(.) %in% c("temp", "prop_mpa")]) %>%
  nest() %>%
  ungroup() %>%
  mutate(placement_id = 1:nrow(.)) %>%
  unnest(cols = data)

write_rds(placement_experiments, file = file.path(results_dir, "placement_experiments.rds"))

write_rds(state_experiments, file = file.path(results_dir, "state_experiments.rds"))

tmp <-
  map(state_experiments$starting_conditions,
      ~ map_df(.x,  ~ sum(.x[[1]]$ssb_p_a) / .x[[1]]$ssb0))

check <-
  tibble(state_id = state_experiments$state_id, tmp = (tmp)) %>%
  unnest(cols = tmp) %>%
  pivot_longer(cols = -state_id,
               names_to = "critter",
               values_to = "depletion")

check %>%
  ggplot(aes(depletion)) +
  geom_histogram() +
  facet_wrap(~ critter)

# run simulations ---------------------------------------------------------
  future::plan(future::multisession, workers = experiment_workers)
  
  experiment_results <-
    vector(mode = "list", length = nrow(placement_experiments))
  
  pb <- progress_bar$new(
    format = "  Running Experiments [:bar] :percent eta: :eta",
    total = nrow(placement_experiments),
    clear = FALSE,
    width = 60
  )
  
  pb$tick(0)
  
  for (p in 1:nrow(placement_experiments)) {
    # memory problem trying to do it all at once so breaking it up a bit
    
    
    # a <- Sys.time()
    tmp <- state_experiments %>%
      ungroup() %>%
      mutate(
        results = future_pmap(
          list(
            starting_conditions = starting_conditions,
            proc_starting_conditions = proc_starting_conditions,
            fauna = fauna,
            fleets = fleet
          ),
          run_mpa_experiment,
          placement_strategy = placement_experiments$placement_strategy[p],
          prop_mpa = placement_experiments$prop_mpa[p],
          critters_considered = placement_experiments$critters_considered[p],
          placement_error = placement_experiments$placement_error[p],
          resolution = resolution,
          patch_area = patch_area,
          .progress = FALSE,
          .options = furrr_options(seed = 42)
        )
      )
    
    
    
    tmp$results <-
      purrr::set_names(tmp$results, state_experiments$state_id)
    # Sys.time() - a
    
    experiment_results[[p]] <- tmp$results
    
    pb$tick()
    
    
  } # close p loop
  
  future::plan(future::sequential)
  
  write_rds(experiment_results,
            file = file.path(results_dir, "experiment_results.rds"))
  
  
} else {
  experiment_results <-
    read_rds(file = file.path(results_dir, "experiment_results.rds"))
  
placement_experiments <- read_rds(file = file.path(results_dir, "placement_experiments.rds"))

state_experiments <- read_rds(file = file.path(results_dir, "state_experiments.rds"))

  
}

# process simulations -----------------------------------------------------


results <- placement_experiments %>%
  mutate(temp = experiment_results,
         state_id = map(experiment_results, names)) %>%
  unnest(cols = c(temp, state_id)) |>
  mutate(fauna = map(temp, c("results", "fauna")),
         fleets = map(temp, c("results", "fleets"))) |>
  select(-temp) |>
  mutate(id = glue::glue("placement_id-{placement_id}_state_id-{state_id}")) |> 
  select(id, everything())


fauna_results <- results |>
  select(-fleets) |>
  unnest(cols = fauna)

fauna_benchmark <- fauna_results |> 
  filter(prop_mpa == 0)

fleet_results <- results |>
  select(-fauna) |>
  unnest(cols = fleets)

fleet_benchmark <- fleet_results |> 
  filter(prop_mpa == 0)

benchmark <- results |> 
  filter(prop_mpa == 0) |> 
  group_by(id) |> 
  nest(.key="control")

# create results ----------------------------------------------------------

# first, calculate no-MPA baseline for metrics in question

tmp <- results |>
  filter(prop_mpa > 0) |>
  group_by(id, placement_id, state_id, prop_mpa) |>
  nest(.key = "treatment") |>
  arrange(id) |>
  left_join(benchmark, by = "id")

calculate_outcomes <- function(treatment, control){
  
  
  # calculate treatment
  
  
  treatment_fauna <- treatment$fauna[[1]] |>
    group_by(critter, step) |>
    summarise(
      biomass = sum(b),
      spawning_biomass = sum(ssb),
      numbers = sum(n),
      mpa_biomass = sum(b[mpa]),
      fished_biomass = sum(b[!mpa])
    ) |> 
    ungroup() |> 
    pivot_longer(-c(critter, step)) |> 
    mutate(fleet = "nature")
  
  
  treatment_fleets <- treatment$fleets[[1]] |>
    group_by(fleet, critter, step) |>
    summarise(
      catch = sum(catch, na.rm = TRUE),
      revenue = sum(revenue, na.rm = TRUE),
      effort = sum(effort, na.rm = TRUE),
      mean_cpue = mean(cpue, na.rm = TRUE)
    ) |> 
    ungroup() |> 
    pivot_longer(-c(fleet, critter, step))
  
  
  # calculate benchmark
  
  mpas <- treatment$fauna[[1]] |> 
    filter(step == max(step)) |> 
    select(x,y, mpa) |> 
    unique()
  
  control_fauna <- control$fauna[[1]] |> 
    select(-mpa) |> 
    left_join(mpas, by = c("x","y")) |> 
    group_by(critter, step) |>
    summarise(
      biomass = sum(b),
      spawning_biomass = sum(ssb),
      numbers = sum(n),
      mpa_biomass = sum(b[mpa]),
      fished_biomass = sum(b[!mpa])
    ) |> 
    ungroup() |> 
    pivot_longer(-c(critter, step), values_to = "control_value") |> 
    mutate(fleet = "nature")
  
  
  control_fleets <- control$fleets[[1]] |>
    select(-mpa) |> 
    left_join(mpas, by = c("x","y")) |> 
    group_by(fleet, critter, step) |>
    summarise(
      catch = sum(catch, na.rm = TRUE),
      revenue = sum(revenue, na.rm = TRUE),
      effort = sum(effort, na.rm = TRUE),
      mean_cpue = mean(cpue, na.rm = TRUE)
    ) |> 
    ungroup() |> 
    pivot_longer(-c(fleet, critter, step), values_to = "control_value")
  
  
  
  treatment_outcomes <- treatment_fauna |> 
    bind_rows(treatment_fleets)
  
  
  control_outcomes <- control_fauna |> 
    bind_rows(control_fleets)
  
  outcomes <- treatment_outcomes |> 
    left_join(control_outcomes, by = c("step", "critter", "fleet","name")) |> 
    mutate(percent_mpa_effect = value / control_value - 1)
  
  # out <- list(base_fauna = base_fauna, base_fleet = base_fleets)
  
  return(outcomes)
}


calculate_gradients <- function(x, prop_mpa){
  
 fauna <-  x$fauna[[1]]
 
 fleets <-  x$fleets[[1]]
 
 
 if (prop_mpa >0 & prop_mpa < 1){
 
   biomass_gradient <- fauna |> 
     group_by(critter, step) |> 
     nest() |> 
     mutate(biomass_gradient = map_dbl(data, ~ lm(log(b) ~ distance_to_mpa_edge, data = .x |> filter(!mpa))$coefficients[[2]])) |> 
     select(-data) |> 
     mutate(fleet = "nature")
   
   cpue_gradient <- fleets |>
     group_by(fleet, critter, step) |>
     nest() |>
     mutate(cpue_gradient = map_dbl(data, ~ lm(log(cpue) ~ distance_to_mpa_edge, data = .x |> filter(!mpa))$coefficients[[2]])) |>
     select(-data)
   
   effort_gradient <- fleets |> 
     group_by(fleet, critter, step) |> 
     nest() |> 
     mutate(effort_gradient = map_dbl(data, ~ lm(log(effort) ~ distance_to_mpa_edge, data = .x |> filter(!mpa))$coefficients[[2]])) |> 
     select(-data) 
   
   fleet_gradients <- cpue_gradient |> 
     left_join(effort_gradient, by = c("step","fleet","critter"))
   
   out <- fleet_gradients |> 
     bind_rows(biomass_gradient) |> 
     ungroup()
 } else {

   out <- data.frame(step = NA, critter = NA, fleet = NA)
 }
 
 
 return(out)
  
}

tmp <- tmp |> 
  ungroup() |> 
  # slice(1) |> 
  mutate(outcomes = map2(treatment,control, calculate_outcomes, .progress = TRUE),
         gradients = map2(treatment,prop_mpa, calculate_gradients, .progress = TRUE))

# then, calculate MPA outcomes


state_variables <- state_experiments |> 
  select(state_id, data) |> 
  unnest(cols = data) |> 
  select(state_id,kiss, mpa_response, starts_with("sigma"), starts_with("spatial"), fleet_model) |> 
  unique() |> 
  mutate(state_id = as.character(state_id))


mpa_outcomes <- tmp |> 
  select(-treatment,-control,-gradients) |> 
  unnest(cols = outcomes) |> 
  ungroup() |> 
  left_join(state_variables, by = "state_id") |> 
  left_join(placement_experiments, by = c("placement_id", "prop_mpa"))

mpa_gradients <- tmp |> 
  select(-treatment,-control,-outcomes) |> 
  unnest(cols = gradients) |> 
  ungroup()

slopes <- mpa_gradients |> 
  group_by(fleet, prop_mpa) |> 
  summarise(real_neg_slope = mean(cpue_gradient < -0.025, na.rm = TRUE),
            real_pos_slope = mean(cpue_gradient > 0.025, na.rm = TRUE))

mpa_outcomes <- mpa_outcomes |> 
  left_join(mpa_gradients,by = join_by(id,placement_id, state_id, prop_mpa, critter, step, fleet))

# make plots --------------------------------------------------------------

fauna_results |> 
  filter(between(prop_mpa, 0.2,0.3), !mpa) |> 
  ggplot(aes(total_mpa_distance, b, color = prop_mpa )) + 
  geom_point() + 
  facet_wrap(~critter) + 
  scale_color_viridis_c()



mpa_outcomes |> 
  filter(name %in% c("biomass", "catch","mean_cpue")) |> 
  ggplot(aes(percent_mpa_effect)) + 
  geom_histogram() +
  facet_wrap(~name)

spillover_vs_catch_plot <- mpa_outcomes |>
  filter(!is.na(cpue_gradient)) |>
  filter(name %in% c("catch")) |>
  ggplot(aes(cpue_gradient, pmin(1.5, percent_mpa_effect), color = prop_mpa)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 3, alpha = 0.33) +
  scale_x_continuous(
    name = "Spillover",
    breaks = c(-3.5, 0, 2),
    labels = c(
      "CPUE higher near MPA",
      'CPUE equal everywhere',
      "CPUE higher far from MPA"
    )
  ) +
  scale_y_continuous(name = "Change in Total Catch", labels = scales::percent) +
  scale_color_viridis_c(
    name = "% area in MPA",
    labels = scales::percent,
    option = "magma",
    guide = guide_colorbar(
      barwidth = unit(12, "lines"),
      frame.colour = "black"
    )
  ) +
  theme(legend.position = "top") +
  coord_cartesian(clip = "off")

spillover_vs_biomass_plot <- mpa_outcomes |>
  filter(!is.na(biomass_gradient)) |>
  filter(name %in% c("biomass")) |>
  ggplot(aes(biomass_gradient, pmin(1.5, percent_mpa_effect), color = prop_mpa)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 3, alpha = 0.33) +
  scale_x_continuous(
    name = "Spillover",
    breaks = c(-1, 0, 0.75),
    labels = c(
      "CPUE higher near MPA",
      'CPUE equal everywhere',
      "CPUE higher far from MPA"
    )
  ) +
  scale_y_continuous(name = "Change in Total Biomass", labels = scales::percent) +
  scale_color_viridis_c(
    name = "% area in MPA",
    labels = scales::percent,
    option = "magma",
    guide = guide_colorbar(
      barwidth = unit(12, "lines"),
      frame.colour = "black"
    )
  ) +
  theme(legend.position = "top") +
  coord_cartesian(clip = "off")


mpa_outcomes |> 
  filter(!is.na(effort_gradient), placement_strategy == "area") |> 
  filter(name %in% c("catch")) |> 
  ggplot(aes(effort_gradient, pmin(2,percent_mpa_effect), color = prop_mpa)) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  geom_point(size = 3, alpha = 0.75) + 
  facet_wrap(~name) + 
  scale_x_continuous(name = "Change in Effort per unit distance from Nearest MPA border", labels = scales::percent) + 
  scale_y_continuous(name = "Change in total catch caused by MPA", labels = scales::percent) + 
  facet_wrap(~critter, scales = "free") + 
  scale_color_viridis_c(name = "% area in MPA",labels = scales::percent)

 mpa_outcomes |> 
  filter(!is.na(biomass_gradient)) |> 
  filter(name %in% c("biomass")) |> 
  ggplot(aes(biomass_gradient, pmin(2,percent_mpa_effect), color = prop_mpa)) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  geom_point(size = 3, alpha = 0.75) + 
  facet_wrap(~name) + 
  scale_x_continuous(name = "Change in CPUE per unit distance from Nearest MPA border", labels = scales::percent) + 
  scale_y_continuous(name = "Change in total biomass cased by MPA", labels = scales::percent) + 
  facet_wrap(~critter, scales = "free") + 
  scale_color_viridis_c(name = "% area in MPA",labels = scales::percent)

 mpa_outcomes |> 
   filter(!is.na(biomass_gradient), placement_strategy == "area") |> 
   filter(name %in% c("mpa_biomass")) |> 
   ggplot(aes(biomass_gradient, pmin(2,percent_mpa_effect), color = prop_mpa)) + 
   geom_hline(yintercept = 0) + 
   geom_vline(xintercept = 0) +
   geom_point(size = 3, alpha = 0.75) + 
   facet_wrap(~name) + 
   scale_x_continuous(name = "Change in CPUE with distance from MPA", labels = scales::percent) + 
   scale_y_continuous(name = "Change in biomass inside MPA caused by MPA", labels = scales::percent) + 
   facet_wrap(~critter, scales = "free") + 
   scale_color_viridis_c(name = "% area in MPA",labels = scales::percent)
 
 plots <- ls()[str_detect(ls(), "_plot$")]
 
 purrr::walk(plots, \(x) ggsave(filename = file.path(plot_dir,glue::glue("{x}.png")), plot = get(x), width = plot_width, height = plot_height))
 