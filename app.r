## Installation of Packages
setwd("Z:/Daniel/MethylBrowsR/MethylBrowsR_2022-06-29")
my_packages <- c("shiny", "shinyWidgets", "shinythemes", "data.table", "shinyalert") 
## Find the non-installed ones 
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
## install them
if(length(not_installed)) install.packages(not_installed)                             

## Load packages 
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(data.table)
library(shinyalert)
library(stringr)
library(ggplot2)

timeoutSeconds <- 300

require(RColorBrewer)
 cols = brewer.pal(12, "Paired")
# cols  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
#                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


source("data.r")
source("webster_plot_code_DLM.r")

# Define UI for MethylBrowsR app ----
ui <- fluidPage(
			
  # App title ----
  headerPanel("MethylBrowsR"),

  tabsetPanel(
  # Tab 1: Age by DNA Methylation
	tabPanel("Age by DNA Methylation", fluid=TRUE, 

  # Sidebar panel for inputs ----
	sidebarPanel(			   
			   sliderInput("range", 
						   label = "Age:", 
						   min = round(min(samps$age)), 
						   max = round(max(samps$age)), 
			   			   value = c(round(min(samps$age)), round(max(samps$age))),
			   			   	step=5),
						   			   
			   checkboxGroupInput("sex", 
							label="Sex:", 
						    choices = c("Male" = "Male",
										"Female" = "Female"), 
				            selected=c("Male", "Female")),
				
				searchInput("cpgsearch",
							label = "CpG ID Search",
							value = NULL,
							placeholder = "e.g. cg05575921",
							btnSearch="Search",
							btnReset="Reset"
							),
								selectInput("scale",
							label = "Scale",
							choices = c("M-values" = "M-values", 
							           "Beta-values" = "Beta-values"),
							selected = "Beta-values"
							)
				),
			 
  # Main panel for displaying outputs ----
	mainPanel(
	plotOutput("plot1")
			 ) # Close mainPanel
			), # Close TabPanel1
	
	
	tabPanel("Regional Plot and Annotation", fluid=TRUE, 

  # Sidebar panel for inputs ----
	sidebarPanel(	
				selectInput("scale2",
							label = "Scale",
							choices = c("M-values" = "M-values", 
							           "Beta-values" = "Beta-values"),
							selected = "Beta-values"
							),
				searchInput("cpgsearch2",
							label = "Search by CpG",
							value = NULL,
							placeholder = "e.g. cg05575921 ",
							btnSearch="Search",
							btnReset="Reset"
							),

					sliderInput("range2", 
						   label = "Age:", 
						   min = round(min(samps$age)), 
						   max = round(max(samps$age)), 
			   			   value = c(round(min(samps$age)), round(max(samps$age))),
			   			   	step=5),
						   			   
			   checkboxGroupInput("sex2", 
							label="Sex:", 
						    choices = c("Male" = "Male",
										"Female" = "Female"), 
				            selected=c("Male", "Female")),


				searchInput("genesearch",
							label = "Search by UCSC Gene ID",
							value = NULL,
							placeholder = "e.g. AHRR",
							btnSearch="Search",
							btnReset="Reset"
							),
			   
			   numericInput("dist", 
							label = "Distance Flanking Target CpG/Gene:", 
						    value=5000
							),
							
			   selectInput("chromosome", 
						   label="Chromosome:", 
						   choices = c("", 1:22, "X", "Y"),
						   selected=""
						   ),
				
			   numericInput("start", 
							label = "Start Position (Hg19):", 
						    value=NULL
							),
			   
			   numericInput("end", 
							label = "End Position (Hg19):", 
						    value=NULL
							)
				), # Close sidebarPanel
				
  # Main panel for displaying outputs ----
	mainPanel(plotOutput("plot2")
		          ) # Close mainPanel
			)	# Close tabPanel 2
		) # Close TabSetPanel
	) # Close UI	  


# Define server logic to plot various variables against mpg ----
# Reactive data frames must be in here
server <- function(input, output) {
  # # Compute the formula text ----
  # # This is in a reactive expression since it is shared by the
  # # output$caption and output$mpgPlot functions
  # formulaText <- reactive({
    # paste(input$variable, "plotted by sex")
  # })

  # # Return the formula text for printing as a caption ----
  # output$caption <- renderText({
    # formulaText()
  # })
  
    plot1_empty <- reactive({ 
    plot(1,0,col="#ffffff00", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    mtext("Please search for a CpG")
  })

        plot1_disclose <- reactive({ 
    plot(1,0,col="#ffffff00", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    mtext("Too few individuals in group")
  })

	# meth = reactive({
	# files = unique(anno$file_ind[which(anno$Name %in% input$cpgsearch)])
	# dnam = list()
	# for(file in files){
	# dnam[[file]] = readRDS(paste0("datafiles/dat_", file, ".rds"))
 #  }
 #  dnam = do.call("rbind", dnam)
 #  return(dnam)
 #  })

 # Make a reactive data frame based on slider/sex filter
 dat_analysis = reactive({
	# if(input$sex == "All"){ 
	# dat1 = samps[which(samps$age >= input$range[1] & 
					  # samps$age <= input$range[2]),]
	# return(dat1)
	# } else {
	ids = samps[which(samps$age >= input$range[1] & 
	                  samps$age <= input$range[2] & 
					  samps$Sex %in% input$sex),]


  files = unique(anno$file_ind[which(anno$Name %in% input$cpgsearch)])
	meth = list()
	for(file in files){
	meth[[file]] = readRDS(paste0("datafiles/dat_", file, ".rds"))
  }
  meth = do.call("rbind", meth)[input$cpgsearch,]

	if(input$scale=="M-values"){
	meth = meth
	} 
	if(input$scale=="Beta-values"){
	meth = M2Beta(meth)
	}
  dat2 = data.frame(CpG=meth[rownames(ids)], age=ids$age, Sex = ids$Sex)
  dat2$Name = input$cpgsearch
  # if(nrow(dat2)<=10){
  # 	dat2 = dat2[0,]
  # }
	return(dat2)
	#	}
	})
	
	
# Subset CpG data frame based on these IDs
 # cpg_df	= reactive({
	# cpg_dat = cpg[which(rownames(cpg) %in% input$cpgsearch_search), which(colnames(cpg) %in% ids())]
	# return(cpg_dat)
	# })
	
  # nas = reactive({
  # na_cpg = ifelse(length(which(is.na(cpg[input$cpgsearch_search,]))>0), 
                       # which(is.na(cpg[input$cpgsearch_search,])), 
					   # 1:ncol(cpg))
  # return(na_cpg)
  # })
  
  # Can make plots from dat (output)
  # Plot 1: Age by CpG plot
  
  # Correlation, P-value and N text
  # all_cor = reactive({
  # paste0("rAll: ", cor.test(plot1_dat()[,"age"], plot1_dat()[,"cpg")$estimate,
         # "; pAll: ", cor.test(plot1_dat()[,"age"], plot1_dat()[,"cpg"])$p.value,
		 # "; NAll: ", nrow(na.omit(plot1_dat()[,c("cpg", "age")])))
					# })

 # male_cor = reactive({
  # paste0("rAll: ", cor.test(plot1_dat()[which(plot1_dat()$sex=="M"), "age"], plot1_dat()[which(plot1_dat()$sex=="M"), "cpg"])$estimate, 
         # "; pAll: ", cor.test(plot1_dat()[which(plot1_dat()$sex=="M"), "age"], plot1_dat[which(plot1_dat()$sex=="M"), "cpg"])$p.value,
		 # "; NAll: ", nrow(na.omit(plot1_dat()[which(plot1_dat()$sex=="M"), c("cpg", "age")])))
					# })					

# female_cor = reactive({
  # paste0("rAll: ", cor.test(plot1_dat()[which(plot1_dat()$sex=="F"), "age"], plot1_dat()[which(plot1_dat()$sex=="F"), "cpg"])$estimate, 
         # "; pAll: ", cor.test(plot1_dat()[which(plot1_dat()$sex=="F"), "age"], plot1_dat[which(plot1_dat()$sex=="F"), "cpg"])$p.value,
		 # "; NAll: ", nrow(na.omit(plot1_dat()[which(plot1_dat()$sex=="F"), c("cpg", "age")])))
					# })					
  plot1 = reactive({
  ggplot(data=dat_analysis(), aes(x=age, y=CpG)) + 
  geom_point() + 
  geom_density_2d(colour=cols[9]) + 
  geom_smooth(method="loess", aes(color=Sex)) + 
  # geom_smooth(method="loess") + 
  scale_color_manual(values=cols[c(4,8)]) + 
  ylab(paste0(input$cpgsearch, " DNA methylation")) + 
  xlab("Age") # + 
  # labs(caption = all_cor())
  # plot(age(), cpg_out(), xlab="Age", ylab=paste0(input$cpgsearch, " methylation"))
  })
  
#   output$plot1 = renderPlot({
#   tryCatch(
#   { plot1() }, # Code that might error goes here, between the { }
#   error = function(e) {""} # Leave this line as-is.
# )
#   })
  
  observeEvent(input$cpgsearch,{
    if (input$cpgsearch == '') {
      output$plot1 <- renderPlot({plot1_empty()})
    }  
    else if (input$cpgsearch %in% dat_analysis()$Name) {
      output$plot1 <- renderPlot({plot1()})
    }
    })

observeEvent(nrow(dat_analysis()),{
    if(nrow(dat_analysis()) <=10){
      output$plot1 <- renderPlot({plot1_disclose()})
    }
    }) 
 

  
  # Plot 2: Regional Plot

	
  # Set chr/start/end if gene option is used
  plot2_chr = reactive({
  if(toupper(input$genesearch) %in% probe.features$gene_name & input$genesearch!="") {
  dat7 = probe.features[which(probe.features$gene_name == toupper(input$genesearch)), "chr"][1]
  }  
  return(dat7)
  })
  
  plot2_start = reactive({
  if(toupper(input$genesearch) %in% probe.features$gene_name & input$genesearch!=""){
  dat8 = min(probe.features[which(probe.features$gene_name == toupper(input$genesearch)), "pos"])
  } 
  return(dat8)
  })
  
  plot2_end = reactive({
  if(toupper(input$genesearch) %in% probe.features$gene_name & input$genesearch!=""){
  dat9 = max(probe.features[which(probe.features$gene_name == toupper(input$genesearch)), "pos"])
  } 
  return(dat9)
  })


plot2_sites = reactive({
	if(input$genesearch =="" & input$cpgsearch2 == "") {
		sites = rownames(probe.features)[which(probe.features$chr == input$chromosome & 
		                                       probe.features$pos >= input$start & 
		                                       probe.features$pos <= input$end)]
		} else 
	if(input$genesearch !="") {
		sites = rownames(probe.features)[which(probe.features$chr == plot2_chr() & 
		                                       probe.features$pos >= plot2_start() & 
		                                       probe.features$pos <= plot2_end())]
		} else 
	if(input$cpgsearch2 !="") {
		chr = probe.features$chr[which(probe.features$Name==input$cpgsearch2)]
		min_dist = probe.features$pos[which(probe.features$Name==input$cpgsearch2)] - input$dist
		max_dist = probe.features$pos[which(probe.features$Name==input$cpgsearch2)] + input$dist
		sites = rownames(probe.features)[which(probe.features$chr == chr & 
		                                       probe.features$pos >= min_dist & 
		                                       probe.features$pos <= max_dist)]
		}
		return(sites) 
	})

# Set dataset for plot2 (beta vs m-vals and subset by age/sex)
  plot2_data = reactive({
  id = samps[which(samps$age >= input$range2[1] & 
	                 samps$age <= input$range2[2] & 
						       samps$Sex %in% input$sex2),]
  
  files = unique(anno$file_ind[which(anno$Name %in% plot2_sites())])
	meth = list()
	for(file in files){
	meth[[file]] = readRDS(paste0("datafiles/dat_", file, ".rds"))
  }
  meth = do.call("rbind", meth)

	if(input$scale2=="M-values"){
	dat6 = meth[,rownames(id)]
	} 
	if(input$scale2=="Beta-values"){
	dat6 = M2Beta(meth[,rownames(id)])
	}
	if(ncol(dat6)>10){
	return(dat6)
  }
	})	
  
  plot2 = reactive({
  if(input$cpgsearch2 == "" & input$genesearch == ""){
  drawRegionalPlot(probe_id = NA,
                chr = input$chromosome,
                start = input$start, #region start 
                end= input$end,  #region end
                title = paste0("chr", input$chromosome, ":", input$start, "-", input$end, " regional plot"),
                sig.cpgs= NA,
								dist=0,
                dataset=plot2_data())
				} 
   if(input$cpgsearch2 %in% rownames(plot2_data())){
	drawRegionalPlot(probe_id = input$cpgsearch2,
                chr = NA,
                start = NA, #gene start 
                end= NA,  #gene end
                title = paste0(input$cpgsearch2, " regional plot"),
                sig.cpgs= NA,
								dist=input$dist,
                dataset=plot2_data())
				} 
	if(toupper(input$genesearch) %in% probe.features$gene_name & input$genesearch!=""){
	drawRegionalPlot(probe_id = NA,
                chr = plot2_chr(),
                start = plot2_start(), #gene start 
                end= plot2_end(),  #gene end
                title = paste0(toupper(input$genesearch), " regional plot"),
                sig.cpgs= NA,
								dist=input$dist,
                dataset=plot2_data())
				} 
  if(input$cpgsearch2 == "" & input$genesearch == "" & input$chromosome==""){
		plot(1,0,col="#ffffff00", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    mtext("Please select a CpG, Gene or Genomic Region")			
    } 
				})
   
 
  output$plot2 = renderPlot({
  plot2()
  }, height = 900, width = 800)
  
}



shinyApp(ui, server)


#  [1] "cg21649660" "cg03577977" "cg23476957" "cg06974684" "cg04097236"
#  [6] "cg23557545" "cg14493349" "cg18722861" "cg18000445" "cg26568284"
# [11] "cg27361590" "cg07901130" "cg02759911" "cg01799681" "cg23642061"
# [16] "cg14153064" "cg13562911" "cg22143569" "cg15557662" "cg05446010"
# [21] "cg16867657" "cg24724428" "cg21572722" "cg16323298" "cg25151806"
# [26] "cg16151520" "cg03295988" "cg21401490" "cg19531247" "cg27046318"
# [31] "cg11011854" "cg19530144" "cg09017996" "cg01142742"