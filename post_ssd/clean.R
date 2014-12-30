### ----- Query US EPA Ecotox database from local server -----------------------
require(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname="ecotox_new")
df <- dbGetQuery(con, "
SELECT 
  CASE
    -- unit conversions to ug/L
   WHEN results.conc1_unit = 'ng/ml'::text THEN results.conc1_mean::numeric
   WHEN results.conc1_unit = 'ug/ml'::text THEN results.conc1_mean::numeric * 1000::numeric
   WHEN results.conc1_unit = 'g/L'::text THEN results.conc1_mean::numeric * 1000000::numeric
   WHEN results.conc1_unit = 'mg/L'::text THEN results.conc1_mean::numeric * 1000::numeric
   WHEN results.conc1_unit = 'ug/L'::text THEN results.conc1_mean::numeric
   WHEN results.conc1_unit = 'ng/L'::text THEN results.conc1_mean::numeric * 0.001::numeric   
   WHEN results.conc1_unit = 'ppm' THEN results.conc1_mean::numeric * 1000::numeric
   WHEN results.conc1_unit = 'ppb' THEN results.conc1_mean::numeric
   WHEN results.conc1_unit = 'ppt' THEN results.conc1_mean::numeric * 0.001::numeric
   WHEN results.conc1_unit = 'AI mg/L' THEN results.conc1_mean::numeric * 1000::numeric
   WHEN results.conc1_unit = 'AI ug/L' THEN results.conc1_mean::numeric
   WHEN results.conc1_unit = 'AI ug/ml' THEN results.conc1_mean::numeric * 1000::numeric
   WHEN results.conc1_unit = 'AI ng/L' THEN results.conc1_mean::numeric * 0.001::numeric
   ELSE NULL::numeric
 END as value,
 species.latin_name as species,
 species.subphylum_div as subphyl, 
 refs.reference_number as ref
FROM 
 coredata.tests
 LEFT JOIN coredata.results ON tests.test_id = results.test_id
 RIGHT JOIN coredata.refs ON tests.reference_number = refs.reference_number
 RIGHT JOIN coredata.chemicals ON tests.test_cas = chemicals.cas_number
 RIGHT JOIN coredata.species ON tests.species_number = species.species_number
WHERE 
 -- only Freshwater, Lab, tests
 tests.organism_habitat = 'Water'
 AND tests.test_location = 'LAB' 
 AND tests.media_type IN ('NR', 'FW') 
 AND results.endpoint IN ('EC50', 'EC50/', 'EC50*', 'EC50*/', 'LC50', 'LC50/', 
 'LC50*', 'LC50*/')
 AND tests.test_cas = 2921882
 AND species.kingdom = 'Animalia'
 AND results.effect = 'MOR'
 AND results.obs_duration_mean = '96' 
 AND results.obs_duration_unit = 'h'
 AND results.conc1_unit IN ('ng/ml', 'ug/ml', 'g/L', 'mg/L', 'ug/L', 'ng/L', 'ppm', 'ppb', 'ppt', 'AI mg/L', 'AI ug/L', 'AI ug/ml', 'AI ng/L')
 AND cast_to_num(results.conc1_mean) NOTNULL
                 ")
dbDisconnect(con)
dbUnloadDriver(drv)

# take geometric mean per species
require(plyr)
df_agg <- ddply(df, .(species, subphyl), summarise,
      val = exp(mean(log(value))),
      n = length(value))

write.table(df_agg, '/home/edisz/Documents/Uni/Projects/blog/post_ssd/ssd_data.csv', sep = ';', row.names = FALSE)
