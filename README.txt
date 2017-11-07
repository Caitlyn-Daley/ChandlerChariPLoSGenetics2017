This folder contains all of the data (measurements in .csv files) and R scripts needed to generate all of the figures and analyses described in the paper. Most of the filenames for the supplementary data & scripts should be relatively self-explanatory. Below are more detailed descriptions of the data files and analysis scripts. 


./Data

    - Combined_Final_Data_05_14.csv:
        This file contains all of the phenotypic data on fly wings associated with the
        main experiment in this paper (allelic series & complementation for the sd and vg
        in the SAM and ORE genetic backgrounds). Each row represents a fly wing. The
        columns indicate the genotypes of both parents, the temperature the fly was reared
        at, the sex of the fly, and the experimental block and replicate. Two measures
        of wing morphology/size are given: the area of the wing (see paper for details
        on how this was measured), and a semiquantitative wing score on a scale from 1 to
        10 (again, see paper for details on how this is quantified). (Figures 2 - 5, S2, S4 - S7)
        
    - Fig6C_S10A_sd_wg_PH3_reformatted_by_CHC.csv
    Antibody (Immuno-fluorescence) staining data for sd mutants for cell proliferation (Figure 6, S10)
       Each row represents a single stained imaginal disc.
       Background = Genetic background of strain
       Mutation = mutant allele
       stain1 = marker protein (WG or CT)
       stain2 = marker for proliferation (PH3 or caspase)
       Count = number of PH3 (or caspase) positive cells
       Wingdisc_Area = area of wing pouch as delimited by WG.
       Stdized_Count = scaled estimate of PH3 positive cells.
       logStdized_Count = log transformed and standardized count.

    - Fig6E_S10B_vg_wg_PH3_reformatted_by_CHC.csv
    	Antibody (Immuno-fluorescence) staining data for vg mutants for cell proliferation (Figure 6, S10)
        see above for columns and row information.

    - Fig6D_sd_margin_length_formatted.csv
    	Immuno-fluorescence for the Wingless protein (WG) to examine the proportion of the wing margin that is intact for sd mutants (Figure 6).
       Each row represents a single stained imaginal disc.
       Background = Genetic background of strain
       Mutation = mutant allele
       Total_MarginLength = the length of the total wing margin 
       Cumulative_SegLength = cumulative length of the WG expressing region (as a measure of margin length)
       Margin_Proportion = Cumulative_SegLength/ Total_MarginLength 

    - Fig6F_vg_margin_length_formatted.csv
        Immuno-fluorescence for the Wingless protein (WG) to examine the proportion of the wing margin that is intact for vg mutants (Figure 6)
       see above for columns and row information.

    - GS_Wings_RawArea_Oct2015.csv
     Wing area measures for the male hemizygous progeny from crosses of the allelic series of sd to 16 DGRP lines. (Figure S3)
     Each row represents a single measured wing.
     allele = the allele of sd used in the experiment
     Line = the DGRP line used as a male parent
     TotalArea = TotalArea of wing
     Semiqt_Raw/Num = semi-quantitative (ordinal) score as a letter or number

    - BackgroundGenotypes.csv
     Genotypic calls from the genotyping chip to examine the extent of introgression of each mutant (Figure S1).

    - JH_2017_vg_temp_background_wing_length_long.csv
        Raw data on fly wing length from the additional temperature manipulation experiment for two vg alleles (vg[1], vg[83b27]) (Figure S8)
    Each row represents a single measured wing.
    Background = wild type genetic background
    Sex = sex
    Genotype = allele of vg (wt, 1, 83b27)
    Vial = replicate vial number
    Temp = rearing temperature (C)
    length = measured wing length


More info below about scripts


/Allelic Series/:

    --Figure2_all_allelic_series_fulldata_area.R
        Generates Figure 2 in the paper (allelic series plots using wing area as the
        phenotypic measure of wing morphology)
        
    --Figure2_all_allelic_series_fulldata_semiqt.R
        Generates Supplementary Figure 2 (allelic series plots using the semiquantitative
        wing score as the phenotypic measure of wing morphology)
        
/Complementation/:

    --Figures3_4_complementation.R
        Generates Figures 3 and 4 (as two separate PDFs)
        
    --FiguresS4_5_full_complementation.R
        Generates Supplementary Figures 4 and 5 (as two separate PDFs)
        
    --FiguresS6_7_reciprocal_effects.R
        Generates Supplementary Figures 6 and 7 (as two separate PDFs)
        
/Fig6_FigS10_related/:

    --Fig6CE_S10AB_wing_pouch_plots_CHC.R
    	Uses Data from Fig6C_S10A_sd_wg_PH3_reformatted_by_CHC.csv and Fig6E_S10B_vg_wg_PH3_reformatted_by_CHC.csv to generate a single pdf with Figures 6C, 6E and Supplementary Figures 10A and 10B 
    	
    --Fig6D_sd_margin_script.R
    	Uses Data from Fig6D_sd_margin_length_formatted.csv to generate a single pdf with Figure 6D
    
    --Fig6F_vg_margin_script.R
    	Uses Data from Fig6F_vg_margin_length_formatted.csv to generate a single pdf with Figure 6F
    
/Genotype Data/:

    --BackgroundGenotypes.csv
        Contains data on the genetic backgrounds of the fly strains used in the main
        experiment. Each row represents a different strain and the columns represent
        different marker loci. The format for the marker loci is as follows:
          -First digit = chromosome (1 = X, 2 = chromosome 2, 3 = chromosome 3)
          -Second digit = chromosome arm (e.g., 21 = 2L, 22 = 2R, 31 = 3L, 32 = 3L; note
           that because the X chromosome is not divided into two arms, all markers on
           chromosome 1 begin with 10)
          -Remaining digits = genomic coordinates
          
    --FigureS1_make_genotype_plot.R
        Generates Supplementary Figure 1
        
/SeverityTest/:

    --Contains data on experiment in which mutant flies were crossed to DGRP lines.
      See the separate README.txt file within this directory for details on that
      experiment.
      
/vg_Background_temp/:

    --Contains data on experiment in which the effects of rearing temperature and
      genetic background were examined in more detail for vestigial mutants. See paper
      for more information on this experiment.
      
    --JH_2017_vg_temp_background_wing_length_long.csv
        Raw data on fly wings from this experiment.
        
    --FigureS8JH_vg_temperature_effect.R
        Generates multiple plots and performs analyses associated with this experiment.
        One of these plots is used in the paper as Supplementary Figure 8 (note: this
        script does not automatically generate a PDF; graphics are drawn to the default
        plotting device in R and must be saved manually).
        