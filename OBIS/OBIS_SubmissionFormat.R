##Example code for preparing three csv files (Occurrence, DNA Derived Data, extendedMeasurementOrFact) for OBIS submission from the GOTeDNA Sample Template
#This code has been prepared for the dataset "DFO Monthly GRDI Survey eDNA metabarcoding (COI) in the Bay of Fundy 2019-2021", which is publicly available in OBIS: https://obis.org/dataset/83040ea6-0e87-4e4e-937d-1b826349c464

#Install packages and Load libraries (remove hashtag to install)

#install.packages("dplyr")
#install.packages("stringr")
#install.packages("purrr")
#install.packages("tibble")
#install.packages("worrms")
#install.packages("openxlsx")
#install.packages("lubridate")
#install.packages("tidyr")
#install.packages("devtools")
#devtools::install_github("iobis/obistools")
#install.packages("Hmisc")
#install.packages("Biostrings")

library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(worrms)
library(openxlsx)
library(lubridate)
library(tidyr)
library(obistools)
library(Hmisc)
library(Biostrings)

#Set working directory to the file location of the completed GOTeDNA Sample Template v.2.2
setwd("C:/Users/HEADK/Desktop/OBIS_prep")                    #edit to the location of your folder

#Pull each metadata and metabarcoding data sheet
GOT <- read.xlsx("GOTeDNA-2_data.xlsx", sheet = "Sample_Metadata")
GOT_metabar <- read.xlsx("GOTeDNA-2_data.xlsx", sheet = "Sample_Metabarcoding data")

#Standardize materialSampleID values between sheets, in case there are slight formatting differences

GOT <- GOT %>%
  mutate(materialSampleID = as.character(materialSampleID),
         materialSampleID = sub("\\.0$", "", materialSampleID))

GOT_metabar <- GOT_metabar %>%
  mutate(materialSampleID = as.character(materialSampleID),
         materialSampleID = sub("\\.0$", "", materialSampleID))

# If multiple datasets are in the template, remove hashtags to filter for select owner contacts, locations, target genes, etc. specific to each data sheet

#GOT <- GOT %>%
#  filter(
#    ownerContact == "anais.lacoursiere@dfo-mpo.gc.ca",
#   decimalLatitude  >= 44.80 & decimalLatitude <= 45.27,   #Filter for coordinates between Bangor and St. John (to exclude PEI data)
#   decimalLongitude >= -68.78 & decimalLongitude <= -66.07
#  )

# Filter the metabarcoding dataset
#GOT_metabar <- GOT_metabar %>%
#  filter(
#    target_gene == "COI",
#  )


#Combine the two data sheets into one by matching columns by the materialSampleID (sample name)
GOT_joined <- GOT %>%
  left_join(
    GOT_metabar,
    by = "materialSampleID",
    relationship = "many-to-many"
  )

#Fix the values in eventDate - This code standardizes dates in case they change format from excel after being opened in R
#NOTE: If dates are formatted as month/day/year prior to running this code, they will come out as NA and won't show up in the GOTeDNA app

GOT_joined <- GOT_joined %>%
  mutate(
    # Always work from a character version
    eventDate_chr = as.character(eventDate),
    
    # Try to interpret as numeric (Excel serial)
    eventDate_num = suppressWarnings(as.numeric(eventDate_chr)),
    
    # Fix 2-digit year dates like 21/8/19 -> 21/8/2019
    eventDate_chr_2y = if_else(
      str_detect(eventDate_chr, "^\\d{1,2}/\\d{1,2}/\\d{2}$"),
      sub("(\\d{1,2}/\\d{1,2}/)(\\d{2})$", "\\120\\2", eventDate_chr),
      NA_character_
    ),
    
    eventDate = case_when(
      # 1) Excel numeric serials (42293 etc.)
      !is.na(eventDate_num) ~ as.Date(eventDate_num, origin = "1899-12-30"),
      
      # 2) dd/mm/yy we just expanded to dd/mm/20yy
      !is.na(eventDate_chr_2y) ~ as.Date(eventDate_chr_2y, format = "%d/%m/%Y"),
      
      # 3) Normal dd/mm/yyyy strings
      TRUE ~ as.Date(eventDate_chr, format = "%d/%m/%Y")
    )
  ) %>%
  select(-eventDate_chr, -eventDate_num, -eventDate_chr_2y)


# Create unique event and occurrence IDs for samples containing coordinates

GOT_joined <- GOT_joined %>%
  # 1) Remove any rows that are controls (controlType not NA/empty)
  filter(is.na(controlType) | controlType == "") %>%
  
  # 2) Drop any rows missing coordinates
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  ) %>%
  
  # 3) Work within each original eventID.x
  group_by(eventID.x) %>%
  mutate(
    row_suffix = row_number(),
    
    # eventID for every remaining row (all have coords now)
    eventID = paste0(eventID.x, "-", row_suffix),
    
    # occurrenceID for every remaining row (present + absent)
    occurrenceID = if_else(
      is.na(occurrenceID),
      paste0("DFO-GRDI-BoF:", materialSampleID, "-", row_suffix),
      paste0(occurrenceID, "-", row_suffix)
    )
  ) %>%
  ungroup() %>%
  select(-row_suffix)


#Use this code if you have values for seq_run_id and lib_id in a separate file (obtain this information from the group that sequenced your data)
#This code will likely need to be edited dependent on how the separate file is formatted, I will keep it commented out for now for simplicity

#map_raw <- read.csv("OBIS information(COI).csv", stringsAsFactors = FALSE)

#map_long <- bind_rows(
#  map_raw %>%
#    select(seq_run_id, lib_id),
#  map_raw %>%
#    select(seq_run_id = seq_run_id.1,
#           lib_id    = lib_id.1)
#) %>%
#  # drop completely empty rows
#  filter(!is.na(lib_id) & lib_id != "")

#map_clean <- map_long %>%
#  # remove ENEG controls if you *don’t* want them mapped
#  filter(!str_detect(lib_id, "ENEG")) %>%
#  mutate(
#    # take digits from start of string
#    materialSampleID = str_extract(lib_id, "^[0-9]+"),
#    materialSampleID = as.character(materialSampleID)   # keep as character
#  ) %>%
#  distinct(materialSampleID, .keep_all = TRUE)

#map_clean <- map_clean %>%
#  # treat blanks as NA
#  mutate(seq_run_id = na_if(seq_run_id, "")) %>%
#  # fill seq_run_id down within the whole table
#  fill(seq_run_id)


#Set up the Occurrence core dataframe
#Edit details of each column name specific to your metadata (see glossary in SOP)
#If including seq_run_id and lib_id from above, remove hashtags at the bottom of this script

occurrence <- GOT_joined %>%
  mutate(
    materialSampleID = as.character(materialSampleID),
    project_contact = ownerContact,
    bibliographicCitation = bibliographicCitation,
    basisOfRecord = basisOfRecord.x,
    locationID = samplingStation,
    materialSampleID = materialSampleID,
    samp_name = materialSampleID,
    technical_rep_id = "1",
    recordedBy = recordedBy.x,
    geo_loc_name = "Canada: Bay of Fundy",
    nameAccordingTo = "WoRMS",
    country = "Canada",
    datasetID = "DFO-GRDI-BoF",
    decimalLatitude = decimalLatitude,
    decimalLongitude = decimalLongitude,
    minimumDepthInMeters = "2",
    maximumDepthInMeters = "3",
    samp_mat_process = "filtration; Following a protocol adapted from Walsh et al. (2009), all samples were peristaltically filtered within an hour of collection through a 0.22 μm Sterivex PVDF pressure filter (Millipore)",
    size_frac = "0.22",
    samp_size = "1",
    samp_size_unit = "L",
    filtrationType = "peristaltic",
    filter_material = "Polyvinylidene difluoride (PVDF)",
    filter_name = "Sterivex",
    samp_store_sol = "Longmire",
    samp_store_temp = "-20",
    nucl_acid_ext_kit = "DNeasy® Blood & Tissue Kit (Qiagen)",
    nucl_acid_ext = "https://doi.org/dx.doi.org/10.17504/protocols.io.bh75j9q6",
    nucl_acid_ext_modify = "modified DNeasy® Blood & Tissue Kit (Qiagen) and QIAshredder (Qiagen) protocol at the Fisheries and Oceans Canada (DFO) Aquatic Biotechnology Lab (ABL) located at the Bedford Institute of Oceanography (Dartmouth, NS, Canada).",
    eventDate = eventDate,
    occurrenceID = occurrenceID,
    occurrenceStatus = case_when(
      organismQuantity > 0 ~ "present",
      organismQuantity == 0 ~ "absent",
      TRUE ~ NA_character_
    ),
    language = "en",
    samp_collect_device = "Niskin bottle",
    materialSampleID = materialSampleID,
    month = month(eventDate),
    bibliographicCitation = bibliographicCitation,
    year = year(eventDate),
    scientificNameAuthorship = ownerContact,
    occurrenceID = occurrenceID,
    samp_name = materialSampleID,
    project_name = "Monthly GRDI survey",
    scientificName = scientificName,
    env_broad_scale = "marine biome (ENVO:00000447)",
    env_local_scale = "water surface (ENVO:01001191)",
    env_medium = "sea water (ENVO:00002149)",
    targetTaxonomicAssay = "invertebrate",
    assay_name = "COI-1 (metabarcoding)",
    assay_type = "metabarcoding",
    target_gene = target_gene,
    seq_meth = "Illumina MiSeq [OBI_0002003]",
    seq_kit = "MiSeq Reagent Kit v3 (Illumina)",
    ref_db = "BOLD database that came with the barque download",
    otu_seq_comp_appr = "[vsearch](https://github.com/torognes/vsearch/releases) v2.14.2+",
    otu_db = "barque v1.7.2 (Mathon et al. (2021))",
    pcr_primer_forward = "GGWACWGGWTGAACWGTWTAYCCYCC",
    pcr_primer_reverse = "ACTTTCGTTCTTGATYRA",
    pcr_primer_name_forward = "mlCOIlintF",
    pcr_primer_name_reverse = "jgHCO2198",
    pcr_primer_reference = "https://doi.org/10.1002/ece3.4213 | https://doi.org/10.1186/1742-9994-10-34 | https://doi.org/10.1111/1755-0998.12138",
    identificationRemarks = identificationRemarks,
    checkls_ver = "1.0.2",
    sampleSizeValue = sampleSizeValue,
    sampleSizeUnit = sampleSizeUnit,
    filter_passive_active_0_1 = "1",
    associatedSequences = "https://www.ncbi.nlm.nih.gov/bioproject/1370813",
    samp_category = "sample",
    samp_collec_device = samp_collect_device,
    project_id = "DFO-GRDI-BoF-COI",
    pcr_0_1 = "1",
    platform = "ILLUMINA",
    instrument = "Illumina HiSeq 1500 [OBI_0003386]", 
    tax_assign_cat = "sequence similarity",
    LClabel = NA_character_
  ) %>%
  filter(
    !is.na(minimumDepthInMeters),
    !is.na(maximumDepthInMeters),
    !is.na(target_gene),
    !is.na(scientificName)
  )
# %>%
#  left_join(
#    map_clean, by = "materialSampleID"
#)

#Look up AphiaID (taxonID) from WoRMS by ScientificName

species_lookup <- occurrence %>%
  distinct(scientificName) %>%
  mutate(
    taxonID = map_int(
      scientificName,
      ~ {
        # Keep NAs as NAs
        if (is.na(.x) || .x == "") return(NA_integer_)
        
        id <- tryCatch(
          wm_name2id(.x, marine_only = FALSE),
          error = function(e) NA
        )
        
        if (is.null(id) || length(id) == 0 || is.na(id[1])) {
          NA_integer_
        } else {
          as.integer(id[1])
        }
      }
    )
  )

# Join taxonID into occurrence, but DO NOT drop rows with no match
occurrence <- occurrence %>%
  select(-any_of(c("taxonID", "scientificNameID",
                   "kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "aphiaID"))) %>%  # we'll re-create these
  left_join(species_lookup, by = "scientificName") %>%
  mutate(
    # optional helper flag: which rows had a WoRMS match?
    worms_match = !is.na(taxonID),
    # keep AphiaID explicitly if you want it separate from taxonID
    aphiaID = taxonID
  )


# Build a taxonomy lookup for each unique non-NA AphiaID

# Helper to pull a rank from wm_classification()
.rank_val <- function(classif, rank_name) {
  if (is.null(classif)) return(NA_character_)
  if (!all(c("rank", "scientificname") %in% names(classif))) return(NA_character_)
  v <- classif$scientificname[classif$rank == rank_name]
  if (length(v) == 0) NA_character_ else as.character(v[1])
}

get_taxonomy_for_id <- function(id) {
  classif <- tryCatch(wm_classification(id), error = function(e) NULL)
  rec     <- tryCatch(wm_record(id),        error = function(e) NULL)
  
  tibble(
    taxonID           = id,
    scientificNameID  = if (!is.null(rec) && "lsid" %in% names(rec))
      rec$lsid else NA_character_,
    kingdom           = .rank_val(classif, "Kingdom"),
    phylum            = .rank_val(classif, "Phylum"),
    class             = .rank_val(classif, "Class"),
    order             = .rank_val(classif, "Order"),
    family            = .rank_val(classif, "Family"),
    genus             = .rank_val(classif, "Genus")
  )
}

unique_ids <- occurrence$taxonID %>%
  unique() %>%
  discard(is.na)

taxonomy_lookup <- map_dfr(unique_ids, get_taxonomy_for_id)


#Join WoRMS taxonomy back to occurrence

occurrence <- occurrence %>%
  select(-any_of(c("kingdom", "phylum", "class",
                   "order", "family", "genus",
                   "scientificNameID"))) %>%
  left_join(taxonomy_lookup, by = "taxonID")


#This code is used to extract information from the otus.database.fasta.gz files (barque) to obtain the OTU/ASVids and sequences for input into the DNA-derived data extension for GBIF/OBIS
#Edit code below to "otu_df" to be specific to the bioinformatic output you receive

setwd("C:/Users/HEADK/Desktop/eDNA_ABL data/APC0222/COI-1_BoF_Run02_20-Sep-22/13_otu_database")

# 1. Read the OTU database FASTA
fa <- readDNAStringSet("COI1_mlCOIintF_jgHCO2198.otus.database.fasta.gz")

# 2. Extract headers and sequences
otu_df <- tibble(
  header   = names(fa),
  sequence = as.character(fa)
) %>%
  # remove leading ">"
  mutate(header = str_remove(header, "^>")) %>%
  # parse header into components
  extract(
    col = header,
    into = c("family", "genus", "species_epithet", "otu_id_num", "read_count"),
    regex = "^([^_]+)_([^_]+)_([^_-]+)-otu-([0-9]+)-([0-9]+)$"
  ) %>%
  mutate(
    otu_id        = paste0("otu_", otu_id_num),
    read_count    = as.integer(read_count),
    scientificName = paste(genus, species_epithet)
  ) %>%
  select(family, genus, species_epithet, scientificName,
         otu_id, read_count, sequence)

otu_df


#DNA derived data extension: Link DNA sequences from otu_df

dna_df <- occurrence %>%
  left_join(
    otu_df %>%
      select(
        scientificName,
        seq_id = otu_id,
        dna_sequence = sequence
      ),
    relationship = "many-to-many"
  )

###########################################
#Run Quality Control checks
#NOTE: problems in quality control checks likely occur from repeating occurrenceIDs resulting from empty values at the bottom of a dataframe (happens if there is not a clean merge between the GOTeDNA "sample metadata" and "sample metabarcoding" sheets)

occur <- occurrence #change name of dataframe in case quality control checks fail and you would like to view the original dataframe for problems (Using the code: View("occurrence") )
dna <- dna_df       #if you do not include the dna sequences using the code above, you can use "dna <- occurrence" instead - seems odd but it will be used for the cleaned_dna object later

#Check that the coordinates are in the correct locations
plot_map(occur)
plot_map_leaflet(occur)

#Check that the taxa names match with WoRMS and merge any corrections
worms <- match_taxa(unique(occur$scientificName)) #select "y" and choose option 1 each time

occur <- merge(occur, worms, by="scientificName")
colnames(occur)[colnames(occur) == "scientificNameID.y"] = "scientificNameID"

#Check that all required fields are present in the occurrence table
#Check the uniqueness of the occurrenceID field (Want to = TRUE)

length(occur$occurrenceID) == length(unique(occur$occurrenceID))

#Can view the dataframe and filter to remove empty values from the organismQuantity column if required (remove hashtags)

#View("occur ")

# occur <- occur %>%
#  filter(
#  !is.na(organismQuantity)
#)


#If quality control is clear, select the term names that are specific to each core and extension file
cleaned_occ <- occur %>%
  select(
    occurrenceID, bibliographicCitation, materialSampleID, eventDate, decimalLatitude, decimalLongitude, scientificName,
    organismQuantity, organismQuantityType, sampleSizeValue, sampleSizeUnit,  associatedSequences,  eventID,
    basisOfRecord, locationID,  recordedBy, country, datasetID, occurrenceStatus, minimumDepthInMeters, maximumDepthInMeters,
    language,  month, year, scientificNameAuthorship, taxonID, scientificNameID, kingdom, phylum, class, order, family, genus
  )

cleaned_dna <- dna %>%
  select(
    project_name, dna_sequence, target_gene, pcr_primer_forward, pcr_primer_reverse, samp_name,
    env_broad_scale, env_local_scale, env_medium, samp_mat_process, size_frac, samp_size,
    samp_size_unit, otu_db, seq_kit, otu_seq_comp_appr, pcr_primer_name_forward, pcr_primer_name_reverse,
    pcr_primer_reference,  occurrenceID
  )

#write csv file and submit to OBIS
write.csv(cleaned_occ, "OBIS_GRDI_BoF_COI_occurrence.csv")
write.csv(cleaned_dna, "OBIS_GRDI_BoF_COI_dnaderiveddata.csv")


#Create a dataframe for the extended Measurement or Fact extension

emof <- dna_df %>%
  mutate(
    occurrenceID,
    measurementType = NA_character_, #This column is where all the unique term names from GOTeDNA and FAIRe go
    measurementValue = NA_character_, #
    measurementUnit = NA_character_,
    measurementTypeID = NA_character_, #This is where the URL link to each of the term names go
    measurementValueID = NA_character_,
    measurementUnitID = NA_character_,
    measurementRemarks = NA_character_
  ) %>%
  select(
    occurrenceID, measurementType, measurementValue, measurementUnit, measurementTypeID, measurementValueID, measurementUnitID, measurementUnitID, measurementRemarks
  )


#Map the URL to the origin of each term name

url_map <- c(
  LClabel                  = "https://github.com/GOTeDNA-OBON",
  ownerContact             = "https://github.com/GOTeDNA-OBON",
  samplingStation          = "https://github.com/GOTeDNA-OBON",
  filtrationType           = "https://github.com/GOTeDNA-OBON",
  totalDNAconc             = "https://github.com/GOTeDNA-OBON",
  unitsDNAconc             = "https://github.com/GOTeDNA-OBON",
  dateFiltration           = "https://github.com/GOTeDNA-OBON",
  timeFiltration           = "https://github.com/GOTeDNA-OBON",
  volumeFiltered           = "https://github.com/GOTeDNA-OBON",
  depthWaterTemp           = "https://github.com/GOTeDNA-OBON",
  occurrenceID             = "https://manual.obis.org/darwin_core.html"
)

#You can remove names below such as "seq_run_id" and "lib_id" if you did not include that information earlier
emof <- dna_df %>%
  select(
    seq_id, samp_category, checkls_ver, assay_name, assay_type, targetTaxonomicAssay,
    geo_loc_name, technical_rep_id, project_contact, seq_run_id, lib_id, project_id,
    pcr_0_1, samp_store_sol, samp_store_temp, platform, instrument, tax_assign_cat,
    LClabel, occurrenceID, nucl_acid_ext_kit, filter_material, seq_run_id, lib_id
  ) %>%
  distinct() %>%
  mutate(across(-occurrenceID, as.character)) %>%
  pivot_longer(
    cols = -occurrenceID,
    names_to = "measurementType",
    values_to = "measurementValue"
  ) %>%
  mutate(
    measurementTypeID = dplyr::recode(
      measurementType,
      !!!url_map,
      .default = "https://github.com/FAIR-eDNA/FAIRe_checklist"   #Default URL is FAIRe
    ),
    measurementUnit    = NA_character_,
    measurementValueID = NA_character_,
    measurementUnitID  = NA_character_,
    measurementRemarks = NA_character_
  )

#write csv file and submit to OBIS
write.csv(emof, "OBIS_GRDI_BoF_COI_emof.csv")




