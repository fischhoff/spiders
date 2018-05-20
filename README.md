# Tritrophic interactions between a fungal pathogen, a spider predator, and the blacklegged tick

File organization

The R markdown file “spider.Rmd” completes the analysis and makes the github markdown file spider.md. Raw data files are in folder “raw_data_files”. Working Rdata files, and temporary .csv files, are in folder “Rdata_files”. Tables and figures for the manuscript are in folder “output”. The current draft of the manuscript is “Tritrophic_EE_20180514_proof.pdf”. 

Study objectives

The blacklegged tick Ixodes scapularis is the primary vector for the bacterium causing Lyme disease in eastern North America and for other medically important pathogens. I. scapularis are vulnerable to attack by fungal pathogens and arthropod predators, but the impacts of interactions between biocontrol agents have not been examined. The biocontrol agent Met52®, containing spores of the generalist entomopathogenic fungus Metarhizium brunneum, is used to control blacklegged ticks with efficacy comparable to chemical acaricides. The brush-legged wolf spider Schizocosa ocreata is a predator of I. scapularis. The Tick Project (www.tickproject.org) is the first large-scale use of Met52 against ticks, raising the importance of understanding its effects on non-target species, particularly tick natural enemies. I conducted a field microcosm experiment to assess the compatibility of Met52 and S. ocreata as tick biocontrol agents. The experiment further determined whether Met52 or wolf spiders affected tick behavior as well as tick survival. 

Analysis decisions

The R markdown file includes analysis that informed study design, as well as analysis of data from the experiment. Here I highlight three key decisions in the analysis:

1)	Estimating sample size needed to detect expected effects on tick survival of wolf spiders, Met52, and wolf spider*Met52 interaction. I used a fully crossed design, with four treatments: spider and Met52, spider and water (control for Met52), no spider and Met52, and no spider and water. To estimate how many replicates were needed of each treatment, I conducted a bootstrap power analysis using data from a previous experiment, conducted by a collaborator, involving flat (unfed) nymphs in microcosms with and without wolf spiders. Using the distributions of tick survival in the previous experiment, and assuming certain effects of Met52 and of Met52*spider interaction, the bootstrap randomization estimated the sample size needed to have 80% power to detect the main and interaction effects at the P<0.05 level. Fortunately, the estimated required sample size was achievable, so the experiment could go forward. 

1)	Spider reproductive status: Model comparison was used to test whether spider reproductive status (gravid or without eggs) influenced effects of spiders on tick survival. Based on results of this model comparison, further analysis excluded potential effects of spider reproductive status. 

2)	Evaluating effects of experiment: Model comparison was used to test for effects of treatment on tick survival, spider survival, and tick questing behavior. For each analysis, alternative models were constructed comprising all combinations of predictor variables, including a null (intercept-only) model. Akaike Information Criterion for small sample sizes (AICc) was used to compare the fits of alternative models. If two models differed in their AICc values by a value less than two, the two models were considered to have similar similar levels of support. 


