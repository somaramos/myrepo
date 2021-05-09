A Study of Sex Associated Traits in Brown Rats and Zebrafish
================
Isaac Ramos

-----

**Introduction**: Both data sets are from a large data bank pooled from
a study published in 2018 called *SAGD: a comprehensive sex-associated
gene database from transcriptomes*. The data sets selected were from a
sequencing of sex associated genes derived from the Brown Rat and the
Zebrafish, or *Rattus norvegicus* and *Danio rerio* respectively. Both
data sets come from genes which code for traits in the brain of both
species in order for comparison statistics to hold some weight when
being interpreted. Each dataset contained the variables: “Gene ID”,
which was a long sequence of letters and numbers coordinating to
specific genes that code for the brain tissue, “Symbol” which was a much
shorter way of making distinctions between the genes, “FPKM(M)” and
“FPKM(F)” which stands for Fragments Per Kilobase of transcript per
Million mapped reads (This is a standard measure taken in RNA and DNA
sequencing), “Chromosome” which detailed which chromosome the gene could
be found on, “log2(M/F ratio)” which was the bias each gene was found in
either sex, and “Padj” which is a statistic term meant to reduce Type 1
errors. This data interested me because I love genetics of all kinds. I
like diving into even the most simple of things, such as length of the
DNA read and seeing how it differs from humans, or in this case
Zebrafish. The datasets provide a comprehensive look at the genetics of
sex associated traits, from a data and number standpoint.

-----

``` r
#Keep only the Chromosomes for Zebrafish
FishKaryotype <- c("1","2","3","4","5","6","7","8",
                   "9","10","11","12","13","14","15",
                   "16","17","18","19","20","21","22",
                   "23","24","25","26","27","28",
                   "29","30","X","Y")
#Clean 
Fish_Clean <- Zebrafish %>% filter(Chromosome %in% FishKaryotype) %>%
  drop_na()

#Keep only the Chromosomes for Rats
RatKaryotype <- c("1","2","3","4","5","6","7","8","9","10",
                  "11","12","13","14","15","16","17","18",
                  "19","20","21","22","23","24","25","26",
                  "27","28","29","30","X","Y")
#Clean
Rat_Clean <- Rat %>% filter(Chromosome %in% RatKaryotype) %>%
  drop_na()

#Lets take a look 
head(Fish_Clean)
```

    ## # A tibble: 6 x 9
    ##   `Gene ID` Symbol `SAGD Group` Species `FPKM(M)` `FPKM(F)` Chromosome
    ##   <chr>     <chr>  <chr>        <chr>       <dbl>     <dbl> <chr>     
    ## 1 ENSDARG0… dio2   SAGD_00136   Zebraf…    51.8     18.9    17        
    ## 2 ENSDARG0… gnat1  SAGD_00136   Zebraf…     2.03     7.10   6         
    ## 3 ENSDARG0… galn   SAGD_00136   Zebraf…    22.1     13.1    25        
    ## 4 ENSDARG0… si:ch… SAGD_00136   Zebraf…    35.2     51.9    7         
    ## 5 ENSDARG0… si:ch… SAGD_00136   Zebraf…     0.144    1.12   17        
    ## 6 ENSDARG0… myhb   SAGD_00136   Zebraf…     0.448    0.0452 6         
    ## # … with 2 more variables: `log2(M/F ratio)` <dbl>, Padj <dbl>

``` r
head(Rat_Clean)
```

    ## # A tibble: 6 x 9
    ##   `Gene ID` Symbol `SAGD Group` Species `FPKM(M)` `FPKM(F)` Chromosome
    ##   <chr>     <chr>  <chr>        <chr>       <dbl>     <dbl> <chr>     
    ## 1 ENSRNOG0… Ddx3   SAGD_00113   Rat         26.4    0.0567  Y         
    ## 2 ENSRNOG0… Eif2s… SAGD_00113   Rat         34.0    0.0771  Y         
    ## 3 ENSRNOG0… Kdm5d  SAGD_00113   Rat         11.1    0.0287  Y         
    ## 4 ENSRNOG0… AC239… SAGD_00113   Rat          1.74   0.149   Y         
    ## 5 ENSRNOG0… AC241… SAGD_00113   Rat          3.93   0.00204 Y         
    ## 6 ENSRNOG0… Eif2s3 SAGD_00113   Rat         36.3   54.2     X         
    ## # … with 2 more variables: `log2(M/F ratio)` <dbl>, Padj <dbl>

*These are our 2 data sets. They are from a data base, which examines
sexually selected traits and the deemed significant information behind
them. The way I reshaped the data was by dropping all of the NA values
and filtering out anything not found in the karyotype of both species.*

``` r
#Joining by chromosome
FullData <- full_join(Rat_Clean,
                      Fish_Clean, by = "Chromosome", 
                      suffix = c(".Rat", ".Zebrafish"))

FullData_Clean <- FullData %>% drop_na()
```

*We utilized a full join for the combination of the two data sets,
Rat\_Clean and Fish\_Clean. We chose this method of merging the data
sets because we wanted to look into the differences in the chromosomal
information, per chromosome, of the two species. In order to do this in
the most effective way, we need as much data as possible, thus a full
join was necessary. While the total number of observations deleted by
dropping the NA values was high, proportionally we did not drop many
data points from the data set.*

``` r
FullData_Clean %>% group_by(Chromosome) %>% summarise(Symbol.Zebrafish = n()) %>% arrange(Symbol.Zebrafish)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

    ## # A tibble: 20 x 2
    ##    Chromosome Symbol.Zebrafish
    ##    <chr>                 <int>
    ##  1 18                   289978
    ##  2 11                   319200
    ##  3 15                   391776
    ##  4 12                   392274
    ##  5 17                   392766
    ##  6 14                   414936
    ##  7 19                   422145
    ##  8 20                   430606
    ##  9 13                   459940
    ## 10 16                   497904
    ## 11 9                    514995
    ## 12 6                    736210
    ## 13 8                    893165
    ## 14 10                  1000890
    ## 15 4                   1121225
    ## 16 7                   1138792
    ## 17 5                   1269844
    ## 18 2                   1284843
    ## 19 3                   1367436
    ## 20 1                   2116632

*From the table above, we can determine the amount of different genes
per chromosome in the data for the Zebrafish. The chromosome with the
least amount of sex associated genes is chromosome 18 while the one with
the most is chromosome 1. *

``` r
#Filter for Chromosome 5
Chromosome5 <- FullData_Clean %>%
  filter(Chromosome == "5") %>% 
  select(`log2(M/F ratio).Rat`,`log2(M/F ratio).Zebrafish`, Symbol.Rat,
         Symbol.Zebrafish) %>% 
  arrange(Symbol.Zebrafish) 

#Calculate the mean of the Bias on Sex
Chromosome5 %>% 
  summarize(MeanBias.Fish = mean(`log2(M/F ratio).Zebrafish`),
            MeanBias.Rat = mean(`log2(M/F ratio).Rat`))
```

    ## # A tibble: 1 x 2
    ##   MeanBias.Fish MeanBias.Rat
    ##           <dbl>        <dbl>
    ## 1        0.0319       0.0529

*If the log2 of the male to female ratio results in a negative number,
the trait is more biased towards females. The opposite is true of a
positive value for the log2 of the M/F ratio.*

*Based on the means of both ratios we can see there is a slight bias
towards males, in both Zebrafish and Rats, in the sex-associated genes
for the brains of both species.*

``` r
Statsignificant <- FullData_Clean %>% 
  #Group by Chromosome
  group_by(Chromosome) %>% 
  #Create a variable of the summarized PADJ values
  mutate(Significance.Zebrafish = min_rank(Padj.Zebrafish)) %>%
  #Keep the Species and gene symbol for identification purposes
  select(Symbol.Rat,Symbol.Zebrafish,Chromosome,Significance.Zebrafish) %>%
  arrange(desc(Significance.Zebrafish))

#Take a look at the data 
head(Statsignificant,10)
```

    ## # A tibble: 10 x 4
    ## # Groups:   Chromosome [1]
    ##    Symbol.Rat Symbol.Zebrafish Chromosome Significance.Zebrafish
    ##    <chr>      <chr>            <chr>                       <int>
    ##  1 Tph1       itsn1            1                           25785
    ##  2 Tph1       dnah6            1                           25785
    ##  3 Tph1       mhc1zea          1                           25785
    ##  4 Tph1       atp1a1a.4        1                           25785
    ##  5 Tph1       chmp2bb          1                           25785
    ##  6 Tph1       mcoln1a          1                           25785
    ##  7 Tph1       mief2            1                           25785
    ##  8 Tph1       b4galt1          1                           25785
    ##  9 Tph1       atp1a1a.1        1                           25785
    ## 10 Tph1       si:dkey-28b4.8   1                           25785

*From this data frame we can see the Chromosome which had the largest
Padj values, or the largest p-values were all associated with Chromosome
one of the Zebra fish.*

``` r
Chromosome3 <- FullData_Clean %>% 
#Filter for Chromosome 3
  filter(Chromosome == "3") %>% 
#Take the SD, Variance and MAD of the Bias for Zebrafish
  summarize(SD_MFratio.Zebrafish = sd(`log2(M/F ratio).Zebrafish`),
            Variance.Of.MFRatio.Zebrafish = var(`log2(M/F ratio).Zebrafish`),
            MAD.Zebrafish = mad(`log2(M/F ratio).Zebrafish`))
#Visualized
Chromosome3
```

    ## # A tibble: 1 x 3
    ##   SD_MFratio.Zebrafish Variance.Of.MFRatio.Zebrafish MAD.Zebrafish
    ##                  <dbl>                         <dbl>         <dbl>
    ## 1                0.572                         0.327         0.143

*All three measures of spread showed the data points were all close to
the mean of with the mad being an low value of .1429486. Lets look into
if this is due to the number of data points or are all of the values
actually close to the mean. *

``` r
#Filter for Chromosome 3
FullData_Clean %>% filter(Chromosome == "3") %>% 
#Take the Max and Min of the Ratio 
  summarize(MaxRatio = max(`log2(M/F ratio).Zebrafish`), 
            MinRatio = min(`log2(M/F ratio).Zebrafish`))
```

    ## # A tibble: 1 x 2
    ##   MaxRatio MinRatio
    ##      <dbl>    <dbl>
    ## 1     5.41    -4.95

*Based on the information above, we can see the low standard deviation
comes from a sheer large number of observations, lowering the SD and the
MAD. This can be understood through the Max Bias being 5.41 and the
“Minimal” Bias being -4.96 we say minimal in quotations because the
negative sign represents bias towards females for that locus. *

``` r
FullData_Clean %>% summarize(n_distinct(Symbol.Rat),
                             n_distinct(Symbol.Zebrafish))
```

    ## # A tibble: 1 x 2
    ##   `n_distinct(Symbol.Rat)` `n_distinct(Symbol.Zebrafish)`
    ##                      <int>                          <int>
    ## 1                    17320                          16188

*The *Rattus norvegicus* has more sex associated genes in the brain than
the Zebrafish. We can draw this conclusion due to the Brown Rat having
more unique entries in the Symbol for the gene variable.*

``` r
#
FullData_Clean %>% 
  summarise(Cor_Male_Female_FPKM.Rat = cor(`FPKM(M).Rat`,`FPKM(F).Rat`), 
            Cor_Male_Female_FPKM.Fish = cor(`FPKM(M).Zebrafish`,
                                            `FPKM(F).Zebrafish`))
```

    ## # A tibble: 1 x 2
    ##   Cor_Male_Female_FPKM.Rat Cor_Male_Female_FPKM.Fish
    ##                      <dbl>                     <dbl>
    ## 1                    0.998                     0.997

*Based on the table above we can see there is an immensely strong
correlation between the FPKM (Fragments per Fragments Per Kilobase of
transcript per Million mapped reads) values for both sexes in both
species.*

``` r
#Keep Numeric Variables
FullData_Num <- FullData_Clean %>% select_if(is.numeric)
#Create a custom color
custom_color <- c("#F4EDCA")
#Save the cor as a data frame 
cor(FullData_Num, use = "pairwise.complete.obs") %>% 
  as.data.frame %>% 
# Convert row names to an explicit variable
  rownames_to_column %>% 
  # Pivot so that all correlations appear in the same column
  pivot_longer(-1, names_to = "other_var", 
               values_to = "correlation") %>%
  # Heatmap with geom_tile
  ggplot(aes(rowname, other_var, fill=correlation)) +
  geom_tile() +  
  # Change the scale 
  scale_color_viridis(option = "D") + 
  #Overlay
  geom_text(aes(label = round(correlation,2)), 
            color = custom_color, size = 4,
            check_overlap = TRUE) + 
  #Creat labels 
  labs(title = "Correlation matrix", x = "variable 1", y = "variable 2")
```

![](Presentation_files/figure-gfm/Heatmap-1.png)<!-- -->

*From this Heatmap of the correlation matrix we can see there is very
little correlation between the numeric variables. However there seems to
be a slightly strong negative correlation between Padj.Rat and the ratio
of M/F bias in Rats.*

``` r
#Creat the data for the scatterplot
scat <- FullData_Clean %>%
  #Selected the Species, Chromosome and the Sex Bias ratio 
  select(Species.Rat,Species.Zebrafish,Chromosome,
         `log2(M/F ratio).Rat`,
         `log2(M/F ratio).Zebrafish`)

#Pivot
scat <- scat %>% 
  pivot_longer(cols = c(Species.Rat,Species.Zebrafish), 
               names_to = "remove", values_to = "Species")

#Remove a dummy variable
scat <- scat %>% select(-remove)
#Take a random sample so RStudio can process code
samp1 <- sample_n(scat, 100)

#Create the scatterplot with the sample
scatPlot1 <- samp1 %>% 
  ggplot(aes(`log2(M/F ratio).Rat`,
             `log2(M/F ratio).Zebrafish`)) +
  #Have color coordinate with Species
  geom_point(aes(color = Species)) +
  scale_color_brewer(palette = "Set2") +
  #Creat labels
  labs(title = "A Look into the Bias of Genes in Each Species")


scatPlot1
```

![](Presentation_files/figure-gfm/Sex%20Ratio-1.png)<!-- -->

*A sample of 100 rows of data was taken because if we used the entire
filtered data set the plot would be almost an entire gigabyte. From this
sample we can see there is no large discrepancy between the Male to
Female ratio of sex associated bias between the two species.*

``` r
#Clean up Boxplot Data
BPData <- FullData_Clean %>%
  select(-`Gene ID.Rat`,-`Gene ID.Zebrafish`,
         -Symbol.Rat,-Symbol.Zebrafish,-Padj.Rat,-Padj.Zebrafish,
         -`SAGD Group.Rat`,-`SAGD Group.Zebrafish`)

#Pivot
BPData <- BPData %>% 
  pivot_longer(cols = c(Species.Rat,Species.Zebrafish),
               names_to = "remove", values_to = "Species")

#Perform more cleaning
BPData <- BPData %>%
  select(-remove,-`log2(M/F ratio).Zebrafish`,
         -`log2(M/F ratio).Rat`,-`FPKM(M).Rat`,
         -`FPKM(M).Zebrafish`)

#Take a sample so RStudio does not suffer
samp2 <- sample_n(BPData, 1000)

#Pivot and remove dummy variable
FinalBPData <- samp2 %>%
  pivot_longer(cols = c(`FPKM(F).Rat`,`FPKM(F).Zebrafish`),
               names_to = "D", values_to = "FPKM") %>% 
  select(-D)

#Custom color and labels
BoxP <- FinalBPData %>% 
  ggplot(aes(x = Chromosome, y = FPKM, fill = Species)) +
  geom_boxplot() + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = "FPKM for Females of the Species") + 
  ylim(0,40)


BoxP
```

    ## Warning: Removed 130 rows containing non-finite values (stat_boxplot).

![](Presentation_files/figure-gfm/Boxplot-1.png)<!-- -->

*The multitude of boxplots above show the FPKM values at each chromosome
for each species in the data set. Similar to how we constructed the
above scatterplot, a random sample of the data set was taken in order
for RStudio to actually be able to create the ggplot and keep the plot
at a reasonable size. Despite the few outliers, the majority of the RNA
sequences was below 40 Fragments(per million). The graph looked at
specifically Female individuals of the species. *

``` r
#Take a sample to prevent an Rstudio meltdown 
pca_samp <- sample_n(FullData_Clean, 1000) 


#Clean up the data for a PCA
pcaRNA <- pca_samp  %>% 
  select_if(is.numeric) %>% 
  scale %>% 
  prcomp


names(pcaRNA)
```

    ## [1] "sdev"     "rotation" "center"   "scale"    "x"

``` r
pcaRNA
```

    ## Standard deviations (1, .., p=8):
    ## [1] 1.42636514 1.40882804 1.25315262 1.12388238 0.84854831 0.65142064 0.05106215
    ## [8] 0.01387514
    ## 
    ## Rotation (n x k) = (8 x 8):
    ##                                   PC1         PC2         PC3          PC4
    ## FPKM(M).Rat                0.44612854  0.54330405 -0.07569555  0.008660479
    ## FPKM(F).Rat                0.44599334  0.54360788 -0.07396532  0.008596405
    ## log2(M/F ratio).Rat       -0.06886474 -0.03917958 -0.70004221  0.068958249
    ## Padj.Rat                   0.07925003  0.03385509  0.70069935 -0.031151887
    ## FPKM(M).Zebrafish         -0.53872013  0.44750518  0.02713416 -0.083757726
    ## FPKM(F).Zebrafish         -0.53715965  0.44620535  0.02518616 -0.107981135
    ## log2(M/F ratio).Zebrafish -0.05058647  0.03982036  0.04152974  0.709794698
    ## Padj.Zebrafish             0.09857731 -0.07539377 -0.06830306 -0.686761693
    ##                                    PC5           PC6           PC7
    ## FPKM(M).Rat                0.003429883  0.0008190427 -0.0002319107
    ## FPKM(F).Rat                0.002693572  0.0018877765  0.0003442360
    ## log2(M/F ratio).Rat       -0.004306439 -0.7063165530 -0.0021701723
    ## Padj.Rat                  -0.029318035 -0.7069366132 -0.0020502598
    ## FPKM(M).Zebrafish         -0.031904406 -0.0050004647 -0.7076193951
    ## FPKM(F).Zebrafish         -0.031659688 -0.0098619091  0.7063754026
    ## log2(M/F ratio).Zebrafish -0.699237810  0.0350841244  0.0125507381
    ## Padj.Zebrafish            -0.712845808 -0.0003995106 -0.0119153565
    ##                                     PC8
    ## FPKM(M).Rat                7.070875e-01
    ## FPKM(F).Rat               -7.071243e-01
    ## log2(M/F ratio).Rat       -1.393847e-03
    ## Padj.Rat                   3.385797e-04
    ## FPKM(M).Zebrafish         -9.081741e-05
    ## FPKM(F).Zebrafish          4.797254e-04
    ## log2(M/F ratio).Zebrafish  4.277281e-04
    ## Padj.Zebrafish             2.879153e-04

``` r
#Rotated Data
head(pcaRNA$x)
```

    ##              PC1        PC2         PC3         PC4         PC5        PC6
    ## [1,]  0.28206899 -0.3730228   0.3077228 -0.03572227 -0.13054190  0.1875871
    ## [2,] -1.14823551 -0.9692571 -11.7560614  0.33928199  0.52803661  1.9840923
    ## [3,] -6.39252703  5.6154874   0.4581066 -1.13754434 -0.63629705 -0.2542046
    ## [4,]  0.02572288 -0.1836523   0.2555819 -0.16375486 -0.05759655 -0.2191567
    ## [5,] -0.08961912 -0.5878721  -2.6945241  0.07957696  0.01783805 -0.5506958
    ## [6,] -0.60651858  0.1993260  -0.1234450 -0.24002274 -0.07612138 -0.3684690
    ##               PC7          PC8
    ## [1,] -0.003187944 -0.004067187
    ## [2,]  0.015934264 -0.007763235
    ## [3,] -0.199190487  0.013891870
    ## [4,]  0.005753292 -0.001200586
    ## [5,] -0.006525753 -0.005941069
    ## [6,]  0.033628891 -0.003163800

``` r
pcaRNA_data <- data.frame(pcaRNA$x, 
  Chromosome = pca_samp$Chromosome)

head(pcaRNA_data)
```

    ##           PC1        PC2         PC3         PC4         PC5        PC6
    ## 1  0.28206899 -0.3730228   0.3077228 -0.03572227 -0.13054190  0.1875871
    ## 2 -1.14823551 -0.9692571 -11.7560614  0.33928199  0.52803661  1.9840923
    ## 3 -6.39252703  5.6154874   0.4581066 -1.13754434 -0.63629705 -0.2542046
    ## 4  0.02572288 -0.1836523   0.2555819 -0.16375486 -0.05759655 -0.2191567
    ## 5 -0.08961912 -0.5878721  -2.6945241  0.07957696  0.01783805 -0.5506958
    ## 6 -0.60651858  0.1993260  -0.1234450 -0.24002274 -0.07612138 -0.3684690
    ##            PC7          PC8 Chromosome
    ## 1 -0.003187944 -0.004067187         12
    ## 2  0.015934264 -0.007763235         10
    ## 3 -0.199190487  0.013891870          1
    ## 4  0.005753292 -0.001200586          3
    ## 5 -0.006525753 -0.005941069          5
    ## 6  0.033628891 -0.003163800         11

``` r
#Visualization
ggplot(pcaRNA_data, aes(x = PC1, y = PC2, color = Chromosome)) +
  geom_point() +
  theme(legend.position = "bottom")
```

![](Presentation_files/figure-gfm/PCA-1.png)<!-- -->

``` r
#Rotation Matrix
pcaRNA$rotation
```

    ##                                   PC1         PC2         PC3          PC4
    ## FPKM(M).Rat                0.44612854  0.54330405 -0.07569555  0.008660479
    ## FPKM(F).Rat                0.44599334  0.54360788 -0.07396532  0.008596405
    ## log2(M/F ratio).Rat       -0.06886474 -0.03917958 -0.70004221  0.068958249
    ## Padj.Rat                   0.07925003  0.03385509  0.70069935 -0.031151887
    ## FPKM(M).Zebrafish         -0.53872013  0.44750518  0.02713416 -0.083757726
    ## FPKM(F).Zebrafish         -0.53715965  0.44620535  0.02518616 -0.107981135
    ## log2(M/F ratio).Zebrafish -0.05058647  0.03982036  0.04152974  0.709794698
    ## Padj.Zebrafish             0.09857731 -0.07539377 -0.06830306 -0.686761693
    ##                                    PC5           PC6           PC7
    ## FPKM(M).Rat                0.003429883  0.0008190427 -0.0002319107
    ## FPKM(F).Rat                0.002693572  0.0018877765  0.0003442360
    ## log2(M/F ratio).Rat       -0.004306439 -0.7063165530 -0.0021701723
    ## Padj.Rat                  -0.029318035 -0.7069366132 -0.0020502598
    ## FPKM(M).Zebrafish         -0.031904406 -0.0050004647 -0.7076193951
    ## FPKM(F).Zebrafish         -0.031659688 -0.0098619091  0.7063754026
    ## log2(M/F ratio).Zebrafish -0.699237810  0.0350841244  0.0125507381
    ## Padj.Zebrafish            -0.712845808 -0.0003995106 -0.0119153565
    ##                                     PC8
    ## FPKM(M).Rat                7.070875e-01
    ## FPKM(F).Rat               -7.071243e-01
    ## log2(M/F ratio).Rat       -1.393847e-03
    ## Padj.Rat                   3.385797e-04
    ## FPKM(M).Zebrafish         -9.081741e-05
    ## FPKM(F).Zebrafish          4.797254e-04
    ## log2(M/F ratio).Zebrafish  4.277281e-04
    ## Padj.Zebrafish             2.879153e-04

``` r
rotation_data <- data.frame(pcaRNA$rotation, 
  variable = row.names(pcaRNA$rotation))

#Visualization
arrowRNA <- arrow(length = unit(0.05, "inches"), type = "closed")
color2 <- c("#CC79A7")

ggplot(rotation_data) +
  geom_segment(aes(xend = PC1, yend = PC2),
               x = 0, y = 0,
               arrow = arrowRNA) +
  geom_text(aes(x = PC1, y = PC2, label = variable),
            hjust = 1.25, size = 3, color = color2,
            check_overlap = TRUE) +
  xlim(-1., 1.25) + coord_fixed()
```

![](Presentation_files/figure-gfm/PCA-2.png)<!-- -->

``` r
#Variance explained by each PC from SD 
percent <- 100* (pcaRNA$sdev^2 / sum(pcaRNA$sdev^2))
percent
```

    ## [1] 25.431469064 24.809955403 19.629893575 15.788894977  9.000428001
    ## [6]  5.304360702  0.032591785  0.002406494

``` r
perc_data <- data.frame(percent = percent,
                        PC = 1:length(percent))
#Visualization
ggplot(perc_data, aes(x = PC, y = percent)) + 
   geom_col() + geom_text(aes(label = round(percent, 2)),
                           size = 4, vjust = -0.5) +
  ylim(0,35)
```

![](Presentation_files/figure-gfm/PCA-3.png)<!-- -->

``` r
#A scree plot for a more effective visual
fviz_screeplot(pcaRNA) + 
  geom_text(aes(label = round(percent, 2)), 
            size = 4, vjust = -0.5, check_overlap = TRUE) +
  ylim(0,30)
```

![](Presentation_files/figure-gfm/PCA-4.png)<!-- -->

*About half of the variance can be explained from the first two
principle components. About 99% of the variation comes from six out of
the eight PC’s. The scree plot shows a steady to sharp decreas in
percentage of explained variance moving left to right.*

-----

References:

Meng-Wei Shi, Na-An Zhang, Chuan-Ping Shi, Chun-Jie Liu, Zhi-Hui Luo,
Dan-Yang Wang, An-Yuan Guo, Zhen-Xia Chen, SAGD: a comprehensive
sex-associated gene database from transcriptomes, Nucleic Acids
Research, Volume 47, Issue D1, 08 January 2019, Pages D835–D840,
<https://doi.org/10.1093/nar/gky1040>

-----
