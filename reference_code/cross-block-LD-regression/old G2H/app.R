#install.packages('ggplot2')
library(shiny)
library(shinythemes)
library(shinyjs)
library(knitr)
library(reshape2)
library(plotly)
#library(tidyverse)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DT)
library(shinydashboard)
library(plotly)
library(RColorBrewer)
library(ggpmisc)
library(Rmisc)



#read dataframe----
#Disease, Acronym, System, Prevalence, Heritability
#dis.data <- read.table("XXX")
dis.data<-data.frame(
    stringsAsFactors = FALSE,
    Disease = c("Brain Damage","Disorder of the Vestibular System",
                "Epilepsy Related Disorders",
                "Extrapyramidal Abnormal Movement Disorders",
                "Ocular Musculoskeletal Disorder","Aneurysm","Atherosclerosis",
                "Cardiac Dysrhythmia","Cardiomyopathy",
                "Cerebrovascular Disease","General Hypertension",
                "Non-Rheumatic Heart Disease",
                "Peripheral Vascular Disease","Cardiac Congenital Anomaly",
                "Oro-Facial Congenital Anomaly","Appendiceal Disease",
                "Esophageal Disease",
                "Functional Digestive Disorder","Gastrointestinal Ulcer",
                "Rectal Anal Disorder","Goiter","Thyrotoxicosis",
                "Type II Diabetes Mellitus",
                "Acquired Coagulation Defect","Pernicious or B12 Deficiency Anemia",
                "Biliary Tract Disease","Chronic Liver Disease",
                "Allergic Rhinitis","Asthma",
                "Atopic Contact Dermatitis","Chronic Sinusitis",
                "Crohns Disease","Chronic Tonsillitis Adenoiditis",
                "Cystitis Urethritis","Esophagitis",
                "Eye Inflammation","Fasciitis","Food Allergy",
                "Gastritis Duodenitis","IBS",
                "Inflammatory Spondylopathies","Lupus Erythematosus","Lymphatic Disorder",
                "Osteoarthritis","Osteomyelitis","Pleuritis",
                "Psoriasis Related Disorders",
                "Rheumatic Heart Disease",
                "Rheumatoid Arthritis Related Conditions","Secondary Arthropathy","Thyroiditis",
                "Type I Diabetes Mellitus",
                "Ulcerative Colitis","Upper Respiratory Inflammation",
                "Arthropod-Borne Diseases","Cellulitis",
                "Chronic Upper Respiratory Infection","Ear Infection",
                "Eye Infection","Fungal Infection",
                "Gastrointestinal Infection","General Bacterial Infection",
                "General Viral Infection","Herpes Simplex",
                "Herpes Zoster","Infestation","Influenza",
                "Kidney Infection","Pneumonia","STI",
                "Septicemia","Skin Infection","UTI","Viral Hepatitis B",
                "Viral Warts HPV","Acne","Alopecia",
                "Carbuncle Furuncle","Chronic Pruritus",
                "Dyschromia","Hair Related Disease","Keratosis",
                "Nail Disease","Oral Soft Tissue Disease",
                "Photodermatitis","Rosacea","Seborrheic Dermatitis",
                "Electrolyte Acid-Base Balance Disorder",
                "Gout Related Crystal Arthropathies","Malnutrition",
                "Metabolic Syndrome X","Mixed Hyperlipidemia",
                "Obesity","Overweight","Vitamin Deficiency",
                "Abnormal Spine Curvature",
                "General Osteochondropathy","General Spondylosis Spine Disorder",
                "Hernia","Jaw Disease",
                "Muscle Ligament Disorder","Spinal Stenosis",
                "Synovium Tendon Bursa Disorder",
                "Benign Bone Connective Tissue Neoplasm","Benign Colon Neoplasm",
                "Benign Digestive Neoplasm","Benign Skin Neoplasm",
                "Benign Uterine Neoplasm","Cervical Cancer",
                "Hemangioma Lymphangioma","Lipoma",
                "Non-Melanoma Skin Cancer","ADHD","Adjustment Disorder",
                "Anxiety Phobic Disorder","Autism",
                "Bipolar Disorder","Depression","Eating Disorder","Migraine",
                "Mood Disorder","OCD","PTSD",
                "Personality Disorder","Schizophrenia Related Psychosis",
                "Sleep Disorder","Substance Abuse","Cataract",
                "Conjunctival Disorders","Corneal Disorders",
                "Glaucoma","Optic Neuritis Neuropathy",
                "Refraction Accomodation Disorders",
                "Vitreous Body Disorder","Eustachian Tube Disorder",
                "External Ear Disorders","General Ear Disorder",
                "Hearing Loss","Tympanic Membrane Disorders",
                "Cranial Nerve Disorder",
                "Nerve Root Plexus Disorders","Peripheral Nerve Disorder",
                "Breast Disorder",
                "Disease of the Female Reproductive Organs","Deviated Nasal Septum","Emphysema COPD",
                "Sleep Apnea","Bladder Disorder",
                "Urinary Calculus"),
    Acronym = c("BRAIND",
                "DVS","EPILEPSY","EAM","OMD","ANEURYSM",
                "ATHRSLS","CARDYS","CARDMYO","CEREBVAS",
                "HYPTEN","NONRHHD","PERIVAS","CARCONG","OFCONG",
                "APPEND","ESOPHD","DIGEST","GASTROULC",
                "RECTANAL","GOITER","THYRODTOX","DM2","ACQCOAG",
                "B12DEF","BILLARY","CHRLIVER","ALLGRHIN",
                "ASTHMA","DERM","CHRSINUS","CROHNS",
                "TONSADEN","CYSURETH","ESOPHI","EYEINFM","FASCI",
                "FOODALL","GASDUO","IBS","INFMSPON","LUPUS",
                "LYMPHD","OSTARTH","OSTMYEL","PLEUI",
                "PSORIAS","RHHD","RA","2NDARTHR","THYRODI",
                "DM1","ULCCOL","UPRESP","ARTHROPOD","CELLUI",
                "CHRUPRESP","EARINF","EYEINF","FUNGINF",
                "GASTROINF","BACTINF","VIRINF","HERPSX",
                "HERPZOS","INFESTA","FLU","KIDNEYINF","PNEUMO",
                "STI","SEPSIS","SKININF","UTI","HEPAB",
                "HPV","ACNE","ALOPE","CARBFUR","CHRPRUR",
                "DYSCHROM","HAIRD","KERATO","NAILD","ORALSOFT",
                "PHOTODERM","ROSACEA","SEBDERM","ELECTROBAL",
                "GOUT","MALNUTR","METAX","HYPLIPD",
                "OBESE","OVERWEIGH","VITAMIND","SPINCURV",
                "OSTCHON","SPONSPIN","HERNIA","JAWD","MUSCLIGM",
                "SPINSTE","SYNTEND","BENBONCON","BENCOLON",
                "BENDIG","BENSKIN","BENUTER","CERVCAN",
                "HEMLYM","LIPOMA","NONMELANO","ADHD","ADJM",
                "GAD","ASD","BIP","MDD","EATD","MIGRAINE",
                "MOODD","OCD","PTSD","PERSON","SCHZ",
                "SLEEPD","SUBABUS","CATARACT","CONJUN","CORNEAL",
                "GLAUCOMA","OPTNEURO","REFRACC","VTRBODY",
                "EUSTUBE","EXTEAR","GENEAR","HEARING",
                "TYMPMEM","CRANNERV","NERVROOT","PERINERV",
                "BREAST","FEMALEREP","DEVNASAL","COPD","SLEEPAPN",
                "BLADDER","URINARY"),
    System = c("CNS",
               "CNS","CNS","CNS","CNS","Cardiovascular",
               "Cardiovascular","Cardiovascular","Cardiovascular",
               "Cardiovascular","Cardiovascular",
               "Cardiovascular","Cardiovascular","Development",
               "Development","Development","Development",
               "Digestive","Digestive","Digestive","Endocrine",
               "Endocrine","Endocrine","Hematologic",
               "Hematologic","Hematologic","Metabolic","Immune",
               "Immune","Immune","Immune","Immune","Immune",
               "Immune","Immune","Immune","Immune","Immune",
               "Immune","Immune","Immune","Immune",
               "Immune","Immune","Immune","Immune","Immune",
               "Immune","Immune","Immune","Immune","Immune",
               "Immune","Immune","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Infectious Disease",
               "Infectious Disease","Integumentary","Integumentary",
               "Integumentary","Integumentary",
               "Integumentary","Integumentary","Integumentary",
               "Integumentary","Integumentary","Integumentary",
               "Integumentary","Integumentary","Metabolic",
               "Metabolic","Metabolic","Metabolic","Metabolic",
               "Metabolic","Metabolic","Metabolic",
               "Musculoskeletal","Musculoskeletal","Musculoskeletal",
               "Musculoskeletal","Musculoskeletal",
               "Musculoskeletal","Musculoskeletal",
               "Musculoskeletal","Neoplastic Process","Neoplastic Process",
               "Neoplastic Process","Neoplastic Process",
               "Neoplastic Process","Neoplastic Process",
               "Neoplastic Process","Neoplastic Process",
               "Neoplastic Process","Neuropsychiatric",
               "Neuropsychiatric","Neuropsychiatric","Neuropsychiatric",
               "Neuropsychiatric","Neuropsychiatric",
               "Neuropsychiatric","Neuropsychiatric",
               "Neuropsychiatric","Neuropsychiatric","Neuropsychiatric",
               "Neuropsychiatric","Neuropsychiatric",
               "Neuropsychiatric","Neuropsychiatric",
               "Ophthalmological","Ophthalmological","Ophthalmological",
               "Ophthalmological","Ophthalmological",
               "Ophthalmological","Ophthalmological","Otic","Otic",
               "Otic","Otic","Otic","PNS","PNS","PNS",
               "Reproductive","Reproductive","Respiratory",
               "Respiratory","Respiratory","Urinary","Urinary"),
    Prevalence = c(0.005171,
                   0.0094773,0.0062072,0.0079537,0.0108953,
                   0.0061756,0.0405971,0.0446648,0.0067188,
                   0.0278079,0.1733443,0.0298674,0.0188383,0.0057065,
                   0.0088176,0.0037063,0.0769255,0.0512153,
                   0.0078813,0.0092751,0.0162348,0.0066566,
                   0.0658072,0.0036455,0.0082315,0.0229805,
                   0.0109512,0.1077544,0.062577,0.0948771,0.04732,
                   0.002883,0.0115291,0.0141099,0.0103158,
                   0.0454456,0.0221069,0.0043772,0.033645,0.0125992,
                   0.0047385,0.0025048,0.0070159,0.0675327,
                   0.0070137,0.012282,0.0107591,0.0074803,
                   0.0091382,0.0052582,0.0055944,0.0120637,0.0034584,
                   0.0298561,0.0020516,0.0606542,0.0190802,
                   0.1064849,0.0529861,0.063031,0.0185966,
                   0.0046626,0.055944,0.0114797,0.009969,0.0036336,
                   0.0259272,0.004979,0.0072973,0.0033133,
                   0.0057054,0.0117576,0.0829755,0.0010317,
                   0.0376943,0.0363373,0.0070571,0.0046866,0.0103824,
                   0.0191183,0.012296,0.0575441,0.0191131,
                   0.010608,0.0117629,0.0105946,0.0112052,
                   0.0268936,0.0117959,0.0035433,0.0046909,0.0372233,
                   0.0393757,0.0057076,0.0378309,0.0104257,
                   0.0034885,0.0810053,0.0262594,0.0075516,
                   0.1209356,0.0174436,0.0389287,0.0050975,
                   0.0392641,0.0061757,0.0671609,0.0175143,0.0086675,
                   0.0104041,0.0075297,0.0542616,0.0221253,
                   0.0296068,0.0629991,0.0011044,0.0084859,
                   0.0379399,0.001584,0.0294831,0.0060963,0.0022379,
                   0.0032488,0.0015509,0.0055023,0.007202,
                   0.044735,0.0445796,0.0094011,0.0148412,
                   0.0288803,0.0032548,0.0030796,0.0185929,0.0201871,
                   0.030942,0.0268587,0.0320604,0.0035774,
                   0.0069878,0.0039269,0.0279672,0.0444457,
                   0.1054416,0.0064934,0.0146425,0.0237316,0.0078014,
                   0.0182585),
    h2 = c(0.226,
           0.295,0.541,0.262,0.357,0.695,0.406,0.24,
           0.604,0.217,0.462,0.344,0.253,0.393,0.44,
           0.26,0.292,0.203,0.306,0.078,0.408,0.405,
           0.561,0.642,0.652,0.219,0.363,0.445,0.457,
           0.202,0.523,0.33,0.574,0.313,0.364,0.292,
           0.35,0.599,0.288,0.287,0.646,0.66,0.215,
           0.256,0.22,0.184,0.304,0.554,0.486,
           0.629,0.635,0.531,0.471,0.287,0.639,0.226,
           0.464,0.244,0.2,0.211,0.295,0.409,0.267,
           0.188,0.127,0.463,0.412,0.221,0.391,0.519,
           0.391,0.243,0.227,0.54,0.289,0.501,0.321,
           0.424,0.227,0.43,0.238,0.344,0.325,0.113,
           0.516,0.397,0.372,0.288,0.49,0.539,0.63,
           0.61,0.64,0.576,0.629,0.328,0.491,0.325,
           0.173,0.205,0.268,0.343,0.18,0.218,
           0.173,0.506,0.547,0.162,0.238,0.377,0.038,
           0.52,0.763,0.561,0.432,0.924,0.676,0.579,
           0.569,0.374,0.521,0.657,0.577,0.52,0.562,
           0.38,0.422,0.31,0.346,0.317,0.698,0.593,
           0.611,0.484,0.274,0.272,0.243,0.217,0.298,
           0.249,0.656,0.213,0.166,0.235,0.489,
           0.534,0.474,0.223,0.357),
    Age.of.Onset = c(10L,25L,4L,15L,1L,36L,43L,
                     19L,32L,34L,31L,25L,32L,0L,0L,9L,4L,2L,23L,7L,
                     21L,21L,31L,19L,25L,20L,27L,3L,2L,1L,4L,15L,2L,
                     11L,16L,5L,19L,0L,12L,16L,26L,22L,3L,36L,12L,20L,
                     8L,25L,23L,23L,16L,17L,20L,0L,6L,5L,1L,0L,0L,
                     2L,0L,4L,0L,5L,18L,3L,2L,8L,1L,17L,16L,0L,6L,
                     22L,7L,12L,13L,5L,7L,15L,5L,37L,12L,2L,18L,19L,
                     0L,21L,33L,12L,18L,32L,15L,9L,25L,10L,8L,21L,
                     14L,14L,15L,33L,15L,12L,40L,28L,13L,29L,19L,5L,
                     22L,18L,6L,8L,13L,2L,14L,15L,1L,14L,10L,9L,12L,
                     12L,16L,9L,18L,48L,13L,16L,29L,13L,9L,29L,1L,
                     2L,1L,3L,2L,4L,19L,27L,20L,16L,9L,15L,20L,14L,
                     22L)
)
#Simulate a set of individuals with polygenic scores
set.seed(2021)
num=1000
#variables----
disease.name <- dis.data$Acronym
prev <- dis.data$Prevalence
h2 <- dis.data$h2


#create the pgs dataframe
prs.dat<-rnorm(num*length(disease.name)) %>% matrix(.,nrow=num,ncol=length(disease.name)) %>% data.frame()
colnames(prs.dat)=disease.name




# format disease type choices
category_choices <- dis.data$System %>% unique()
names(category_choices) <- category_choices%>% unique()
category_map <- data.frame(category = dis.data$System,
                           disease = dis.data$Acronym)

#main structure of webpage
ui <- dashboardPage(
    dashboardHeader(title = "My Health Risk Dashboard"),
    dashboardSidebar(
        selectizeInput(inputId ='category_id', 
                       label = "Disease Category", 
                       choices = as.factor(category_choices),
                       options = list(
                           placeholder = 'Select a Category'
                       ),
                       selected='Neuropsychiatric'
                       
        ),
        selectizeInput(inputId ='include', 
                       label = "Include your behavioral and environmental data?", 
                       choices = as.factor(c('Yes','No')),
                       options = list(
                           placeholder = 'Please choose Yes or No'
                       ),
                       selected='No'
                       
        ),
        sliderInput("r2",
                    "Average proportion of genetic variance explained (PGVE) by Polygenic Risk Score",
                    min = 0,
                    max = 1,
                    step = 0.01,
                    value = 0.35),
        numericInput("individual.id", "Individual ID No.", min = 1,max = num, step = 1, value = 100),
        tags$p("Note: \n\r
  1.PGVE is a measure of accuracy PRS on recoverying true genetic risk. It ranges between 0% and 100%.\n\r
  2.Reference data is from the Nature Genetics Paper (Kanix Wang et al.,2017) \n\r
  3.Normal relative genetic risk (RGR) is defined as RGR falling 
           between the first (25%) and the third quartile (75%) of the population; similarly, 5%≤RGR<25% → low risk, RGR<5% → very low risk, 75%≤RGR<95% → high risk, RGR≥95% → very high risk. With a very high relative risk means that the individual is genetically more disposed to the disease than over 95% of all people.")),
    dashboardBody(
        fluidRow(
            valueBoxOutput("traits"),
            valueBoxOutput("high_risk"),
            valueBoxOutput("low_risk")
        ),
        fluidRow(
            box(title = "Relative Genetic Risk Profile [RGR = (your GR/population average GR-1)%]", 
                subtitle='aa',
                solidHeader = T,
                width = 15, collapsible = T,
                plotlyOutput("grisk_plot"))
        ), #row
        fluidRow(
            box(title = "Estimated heritability, onset age and population prevalence for the disease category", solidHeader = T,
                width = 15, collapsible = T,
                plotlyOutput("grisk_plot3"))
        ), #row
        fluidRow(
            box(title = "Genetic Risk Profile [your GR = Probability of getting the disease given your genetics]", solidHeader = T,
                width = 15, collapsible = T,
                plotlyOutput("grisk_plot2"))
        ), #row
        fluidRow(
            box(title = "Risk Components Profile for the traits you have high risk with", solidHeader = T,
                width = 15, collapsible = T,
                plotlyOutput("grisk_plot4"))
        ) #row
    ) # body
    
)


# Define server logic required to draw a histogram
#do it for individual #1000


risk_cols <- brewer.pal(5,name="RdYlGn")

server <- function(input, output) {
    #iid <- input$individual.id
    r2mean<-reactive(input$r2)
    input.category <- reactive(input$category_id)
    id <- reactive(input$individual.id)
    incl <- reactive(input$include)
    rawdat<-reactive({ #crease relative risk df
        #intercept
        #genetic variance explained
        r22mean=r2mean()
        set.seed(1)
        r22=rnorm(length(prev),mean=r22mean,sd=0.05)
        risk.dat=prs.dat
        int=log(prev/(1-prev))
        for(i in 1:length(prev)){
            risk.dat[,i]=sqrt(h2[i]*r22[i]*prev[i]*(1-prev[i]))*prs.dat[,i]+int[i]}
        
        #linear approximation
        #raw.dat[,i]=prev[i]+sqrt(h2[i]*r2[i]*prev[i]*(1-prev[i]))*prs.dat[,i]
        #risk.dat=raw.dat  
        #relative genetic risk
        
        risk.dat=1/(1+exp(-risk.dat))
        raw.dat=risk.dat
        raw.dat})
    grisk_df <- reactive({
        raw.dat<-rawdat()
        raw.dat[,names() %in% category_map$disease[category_map$category == input.category()]]})
    dis_df <- reactive({dis.data[dis.data$System %in% input.category(),]})
    
    
    
    grisk_pct <- reactive({
        raw.dat<-rawdat()
        temp<-raw.dat[,names(raw.dat) %in% category_map$disease[category_map$category == input.category()]]
        grisk_pct <- temp
        iid <- id()
        dis.tmp <- dis_df()
        for(i in 1:ncol(grisk_pct)) {
            grisk_pct[,i] <- (grisk_pct[,i]-dis.tmp$Prevalence[i])/dis.tmp$Prevalence[i]
        }
        
        q5=apply(grisk_pct,2,function(x)quantile(x,0.05))
        q25=apply(grisk_pct,2,function(x)quantile(x,0.25))
        q50=apply(grisk_pct,2,function(x)quantile(x,0.50))
        q75=apply(grisk_pct,2,function(x)quantile(x,0.75))
        q95=apply(grisk_pct,2,function(x)quantile(x,0.95))
        
        risk<-as.numeric(temp[iid,])
        
        grisk_pct <- grisk_pct[iid,] %>% pivot_longer(cols = everything(), names_to = "disease", values_to = "grisk.pct")
        
        grisk_pct$q5 <- q5
        grisk_pct$q25 <- q25
        grisk_pct$q50 <- q50
        grisk_pct$q75 <- q75
        grisk_pct$q95 <- q95
        grisk_pct$risk<-risk
        
        grisk_pct$risk_cat <- NA
        grisk_pct$risk_cat[grisk_pct$grisk.pct <= grisk_pct$q5] <- "Very Low Risk"
        grisk_pct$risk_cat[grisk_pct$grisk.pct > grisk_pct$q5 & grisk_pct$grisk.pct <= grisk_pct$q25] <- "Low Risk"
        grisk_pct$risk_cat[grisk_pct$grisk.pct > grisk_pct$q25 & grisk_pct$grisk.pct <= grisk_pct$q75] <- "Normal Risk"
        grisk_pct$risk_cat[grisk_pct$grisk.pct > grisk_pct$q75 & grisk_pct$grisk.pct <= grisk_pct$q95] <- "High Risk"
        grisk_pct$risk_cat[grisk_pct$grisk.pct >= grisk_pct$q95] <- "Very High Risk"
        
        #grisk_pct <- grisk_pct %>% mutate(risk_cat = case_when(grisk.pct <= q5 ~ "Very Low Risk",
        #                                                       grisk.pct > q5 & grisk.pct <= q25 ~ "Low Risk",
        #                                                       grisk.pct > q25 & grisk.pct <= q75 ~ "Normal Risk",
        #                                                       grisk.pct > q75 & grisk.pct <= q95 ~ "High Risk",
        #                                                       grisk.pct >= q95~ "Very High Risk"))
        grisk_pct$risk_cat <- factor(grisk_pct$risk_cat, levels = c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk"))
        grisk_pct$disease <- factor(grisk_pct$disease)
        grisk_pct
    })
    
    
    output$grisk_plot <- renderPlotly({
        iid<-id()
        dis.tmp <- dis_df()
        grisk_tmp <- grisk_pct()
        dis.tmp <-dis.tmp[order(grisk_tmp$grisk.pct),]
        grisk_tmp <-grisk_tmp[order(grisk_tmp$grisk.pct),]
        grisk_tmp$disease<-dis.tmp$Disease
        #b=c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk")
        grisk_tmp$disease<-factor(grisk_tmp$disease,levels=grisk_tmp$disease)
        grisk_tmp$relative_genetic_risk <- grisk_tmp$grisk.pct*100
        grisk_tmp$risk_category<-factor(grisk_tmp$risk_cat,levels=c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk"))
        grisk_tmp$relative_genetic_risk<-round(grisk_tmp$relative_genetic_risk,4)
        plot_grisk <- ggplot(grisk_tmp) + 
            geom_bar(aes(x = disease, y = relative_genetic_risk, fill = risk_category),stat='identity') +
            labs(x = "Disease", y = "Relative Genetic Risk (compared to others)%", fill = "Relative Risk Group") + 
            theme_bw() + 
            scale_x_discrete(labels= dis.tmp$Acronym,drop=F)+
            scale_fill_manual(values = risk_cols[5:1],
                              labels =c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk"),
                              drop=F)
        ggplotly(plot_grisk)
    })
    
    output$grisk_plot2 <- renderPlotly({
        iid<-id()
        dis.tmp <- dis_df()
        grisk_tmp <- grisk_pct()
        #color=risk_cols[5:1]
        #names(color)=c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk")
        dis.tmp <-dis.tmp[order(grisk_tmp$risk),]
        grisk_tmp <-grisk_tmp[order(grisk_tmp$risk),]
        grisk_tmp$disease<-dis.tmp$Disease
        grisk_tmp$disease<-factor(grisk_tmp$disease,levels=grisk_tmp$disease)
        grisk_tmp$genetic_risk <- grisk_tmp$risk*100
        grisk_tmp$risk_category<-grisk_tmp$risk_cat
        grisk_tmp$genetic_risk<-round(grisk_tmp$genetic_risk,4)
        plot_grisk <- ggplot(grisk_tmp) + 
            geom_col(aes(x = disease, y = genetic_risk, fill = risk_category)) +
            labs(x = "Disease", y = "Genetic Risk (probability)%", fill = "Relative Risk Group") + 
            theme_bw() + 
            scale_x_discrete(labels= as.factor(dis.tmp$Acronym),drop=FALSE)+
            #ylim(0, 100) +
            #geom_hline(yintercept = 50) +
            scale_fill_manual(values = risk_cols[5:1], 
                              labels =c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk"),
                              drop=F
            )
        ggplotly(plot_grisk)
    })
    
    
    output$grisk_plot3 <- renderPlotly({
        
        dis.tmp <- dis_df()
        grisk_tmp <- grisk_pct()
        dis.tmp <-dis.tmp[order(grisk_tmp$grisk.pct),]
        dis.tmp2<-dis.tmp[,c('Acronym','h2','Prevalence','Age.of.Onset')]
        colnames(dis.tmp2)<-c('Disease','Heritability','Prevalence','Onset Age')
        #dis.tmp2$Heritability<-dis.tmp2$Heritability %>% round(.,3)
        #dis.tmp2$Prevalence<-dis.tmp2$Prevalence %>% round(.,4)
        #dis.tmp2$`Onset Age`<-dis.tmp2$`Onset Age` %>% round(.,0)
        dis.tmp2$Disease<-factor(dis.tmp2$Disease,levels=dis.tmp2$Disease)
        dis.tmp2<-gather(dis.tmp2,'property','value',-Disease)
        dis.tmp2$value<-dis.tmp2$value %>% round(.,4)
        plot_grisk <- ggplot(dis.tmp2,aes(x = Disease, y = value, fill = property )) + 
            geom_bar(stat='identity',position = 'dodge') +
            facet_wrap(.~property,scales = 'free_y',nrow=3)
        labs(x = "Disease", y = "Value", fill = NULL) + 
            geom_text(label= dis.tmp$Disease)+
            #theme(axis.text.x = theme_text(angle=0, size=8, hjust=1))
            theme(strip.text.x = element_blank())+
            theme_bw() 
        
        ggplotly(plot_grisk)
    })
    
    
    output$grisk_plot4 <- renderPlotly({
        include.id=incl()
        iid<-id()
        dis.tmp <- dis_df()
        grisk_tmp <- grisk_pct()
        #color=risk_cols[5:1]
        #names(color)=c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk")
        dis.tmp <-dis.tmp[order(grisk_tmp$risk),]
        grisk_tmp <-grisk_tmp[order(grisk_tmp$risk),]
        grisk_tmp$disease<-dis.tmp$Disease
        grisk_tmp$disease<-factor(grisk_tmp$disease,levels=grisk_tmp$disease)
        grisk_tmp$genetic_risk <- grisk_tmp$risk*100
        grisk_tmp$risk_category<-grisk_tmp$risk_cat
        grisk_tmp$genetic_risk<-round(grisk_tmp$genetic_risk,4)
        df=grisk_tmp[c(grisk_tmp$risk_cat == "High Risk"|grisk_tmp$risk_cat == "Very High Risk"),]
        df2=dis.tmp[c(grisk_tmp$risk_cat == "High Risk"|grisk_tmp$risk_cat == "Very High Risk"),]
        
        # Create Data
        if(include.id=='No'){
            df3 <- df[,c('disease','genetic_risk')]
            df3$environmental_risk=(1-sqrt(df2$h2))/sqrt(df2$h2)*df3$genetic_risk
            
        }
        if(include.id=='Yes'){
            df3 <- df[,c('disease','genetic_risk')]
            datatype=c('smoking behavior','drinking behavior','sleeping pattern','nutrition','diet','height','BMI','Metabolic')
            prop=rnorm(length(datatype),1,0.5)
            prop=prop/sum(prop)
            a=ncol(df3)
            for(i in 1:length(datatype)){
                df3[,a+i]=(1-sqrt(df2$h2))/sqrt(df2$h2)*df3$genetic_risk*prop[i]
                colnames(df3)[ncol(df3)]=datatype[i]
            }
            #  df3$environmental_risk=(1-sqrt(df2$h2))/sqrt(df2$h2)*df3$genetic_risk
            
        }
        
        df4<-gather(df3,key='Risk_Component',value='Risk',-disease)
        plot_grisk <- ggplot(df4) + 
            geom_col(aes(x = disease, y = Risk, fill = Risk_Component),stat='identity',position = "stack") +
            labs(x = "Disease", y = "Total Risk (probability)%", fill = "Risk Component Group") + 
            theme_bw() 
        #scale_x_discrete(labels= as.factor(dis.tmp$Acronym),drop=FALSE)+
        #ylim(0, 100) +
        #geom_hline(yintercept = 50) +
        # scale_fill_manual(values = risk_cols[5:1], 
        #                   labels =c("Very Low Risk","Low Risk", "Normal Risk", "High Risk","Very High Risk"),
        #                   drop=F
        # )
        ggplotly(plot_grisk)
    })
    
    
    # value boxes
    output$traits <- renderValueBox({
        aa<- dis_df()
        valueBox(paste0(nrow(aa), " Traits"), 
                 "Traits calculated for this category", icon = icon("file-medical"), color = "blue")
    })
    output$high_risk <- renderValueBox({
        grisk_tmp <- grisk_pct()
        valueBox(paste0(sum(grisk_tmp$risk_cat == "High Risk")+sum(grisk_tmp$risk_cat == "Very High Risk"), " High Risk"), 
                 "No. of traits with (very) high relative genetic risk", icon = icon("exclamation-triangle"), color = "red")
    })
    
    output$low_risk <- renderValueBox({
        grisk_tmp <- grisk_pct()
        valueBox(paste0(sum(grisk_tmp$risk_cat == "Low Risk")+sum(grisk_tmp$risk_cat == "Very Low Risk"), " Low Risk"), 
                 "No. of traits with (very) low relative genetic risk", icon = icon("check-circle-o"), color = "green")
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
