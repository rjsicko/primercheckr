library(openssl)

library(BiocManager)
options(repos = BiocManager::repositories())

BiocManager::install("GenomeInfoDbData")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("MafDb.gnomAD.r2.1.hs37d5")
BiocManager::install("MafH5.gnomAD.v3.1.2.GRCh38")

library(shiny)
library(DT)
library("dplyr")
library(GenomeInfoDb)
library(BSgenome)
library("MafDb.gnomAD.r2.1.hs37d5")
library("MafH5.gnomAD.v3.1.2.GRCh38")
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
mafdb_hg19 <- MafDb.gnomAD.r2.1.hs37d5
mafdb_hg38 <- MafH5.gnomAD.v3.1.2.GRCh38
hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg38 <- BSgenome.Hsapiens.UCSC.hg38

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            tags$b("primerCheck"),
            tags$h4("v0.0.1-alpha1"),
            tags$hr(),
            fileInput("file1", "Load a primer file (csv)", accept = ".csv"),
            #actionButton("exampleData", label = "Load Example Data"),
            width = 4,
            tags$b("To Do"),
            tags$hr(),
            tags$ol(
              tags$li("Help/example page"),
              tags$li("MAF cutoff input to filter tables"), #https://rstudio.github.io/DT/options.html
              tags$li("gnomAD links in the SNP tables"),
              tags$li("Additional population databases?"),
              tags$li("Handle input of single primer"),#maybe blat for coordinates
              tags$li("PCR products tab w/ expected amplicons"),
              tags$li("UCSC browser links for regions"),
              tags$li("Compressed table, where max MAF variant displayed and all others 'under' that. e.g. https://rstudio.github.io/DT/002-rowdetails.html"),
              tags$li("Highlight base in primer row of data corresponds to"),
              tags$li("accept M13 tagged primers (trim tags internally)")) 
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Primers(input)", DT::dataTableOutput("primers")),
                tabPanel("SNPs_gnomAD.v3.1.2(hg38)", DT::dataTableOutput("hg38")),
                tabPanel("SNPs_gnomADv2.1(hg19)", DT::dataTableOutput("hg19"))
            )
        )
    )
)




server <- function(input, output){
  
  
  primer_file <- reactive(
  {
    #observeEvent(input$exampleData,
    #    primer_file <- data.frame (ID  = c("Pompe-E2", "IDUA_E2", "MCADD-E10-F"),
    #                               chr = c("chr17", "chr4", "chr1"),
    #                               FP = c("GCCACCTTACCCCACCTG", "ATGGGACCTCCCCACATTCC", "CCACCAGTTTCTTGTTGCT"),
    #                               RP = c("TGCGTGGGTGTCGATGT", "GAACAAGGGGTCTTCCGAGC", "CCACATTATCTTCTTGCTATCAAAGT")))
    #
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    primer_file <- as.data.frame(read.csv(file$datapath, header = TRUE))
    return(primer_file)
    
  })
 #primer_file <- observeEvent(input$exampleData,{
 #  primer_file <- data.frame (ID  = c("Pompe-E2", "IDUA_E2", "MCADD-E10-F"),
 #                             chr = c("chr17", "chr4", "chr1"),
 #                             FP = c("GCCACCTTACCCCACCTG", "ATGGGACCTCCCCACATTCC", "CCACCAGTTTCTTGTTGCT"),
 #                             RP = c("TGCGTGGGTGTCGATGT", "GAACAAGGGGTCTTCCGAGC", "CCACATTATCTTCTTGCTATCAAAGT"))
 #  return(primer_file)
 #  
 #})
  output$primers <- DT::renderDT(primer_file(),
                                 escape=F,
                                 filter = list(position = 'bottom', clear = FALSE, plain = TRUE),
                                 rownames= F,
                                 extensions = list("ColReorder" = NULL,
                                                   "Buttons" = NULL,
                                                   "FixedColumns" = list(leftColumns=1)),
                                 options = list(
                                   dom = 'BRrltpi',
                                   autoWidth=TRUE,
                                   lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                                   ColReorder = TRUE,
                                   buttons =
                                     list(
                                       'copy',
                                       'print',
                                       list(
                                         extend = 'collection',
                                         buttons = c('csv', 'excel', 'pdf'),
                                         text = 'Download'
                                       ),
                                       I('colvis')
                                     )
                                 )
       )
                                                                  
  
  my_data <- reactive(
  {
     primers_hg38 <- data.frame(primer_id=character(),
                                        seq=character(),
                                        chr=character(),
                                        hg38_pos=integer(),
                                        AF_allpopmax_hg38=integer(),
                                        stringsAsFactors=FALSE)
  
     primers_hg19 <- data.frame(primer_id=character(),
                                        seq=character(),
                                        chr=character(),
                                        hg19_pos=integer(),
                                        AF_afr=integer(),AF_amr=integer(), AF_asj=integer(), AF_eas=integer(), AF_fin=integer(), AF_nfe=integer(), AF_oth=integer(),
                                        stringsAsFactors=FALSE)
             
    progress <- shiny::Progress$new()
       # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 0)
               
    lastrow <- nrow(primer_file())
    firstrow=1
    for (no in (firstrow:lastrow))
    {
       row = primer_file()[no,]
       temp_chr <- row$chr
       temp_FP <- row$FP
       temp_RP <- row$RP
       
       progress$inc(1/lastrow, detail = paste("line ", no))
         
         ###################lets do hg38 first######################
         
         subject_hg38 <- hg38[[temp_chr]]
         products_hg38 <- matchProbePair(temp_FP,temp_RP,subject_hg38)
         amp_start_hg38 = (start(products_hg38))
         fp_end_hg38 = (start(products_hg38) + as.integer(nchar(row$FP)) - 1)
         rp_start_hg38 = ( (start(products_hg38) + width(products_hg38)) - as.integer(nchar(row$RP)) )
         amp_end_hg38 = (start(products_hg38) + width(products_hg38) - 1)
         
         fp_range_hg38 <- GRanges(seqnames=temp_chr, IRanges(start=amp_start_hg38:fp_end_hg38, width=1))
         fp_scores_hg38 <- gscores(mafdb_hg38,fp_range_hg38,pop="AF_allpopmax")
         fp_scores_hg38 <- as.data.frame(fp_scores_hg38)
         
         rp_range_hg38 <- GRanges(seqnames=temp_chr, IRanges(start=rp_start_hg38:amp_end_hg38, width=1))
         rp_scores_hg38 <- gscores(mafdb_hg38,rp_range_hg38,pop="AF_allpopmax")
         rp_scores_hg38 <- as.data.frame(rp_scores_hg38)
         
         #primers_hg38 <- data.frame(primer_id=character(),
         #			seq=character(),
         #			chr=character(),
         #			hg38_pos=integer(),
         #			AF_allpopmax_hg38=integer(),
         #			stringsAsFactors=FALSE)
         primers_hg38 <- bind_rows(primers_hg38,setNames(data.frame(primer_id <- paste0(row$ID,"_F"),
                                                                    seq <- temp_FP,
                                                                    chr <- fp_scores_hg38$seqnames,
                                                                    hg38_pos <- fp_scores_hg38$start,
                                                                    AF_allpopmax_hg38 <- fp_scores_hg38$AF_allpopmax, stringsAsFactors=FALSE),c("primer_id", "seq", "chr","hg38_pos","AF_allpopmax_hg38")))
         #names(primers_hg38) <- c("primer_id","seq","chr","hg38_pos","AF_allpopmax_hg38")
         primers_hg38 <- bind_rows(primers_hg38,setNames(data.frame(primer_id <- paste0(row$ID,"_R"),
                                                                    seq <- temp_RP,
                                                                    chr <- rp_scores_hg38$seqnames,
                                                                    hg38_pos <- rp_scores_hg38$start,
                                                                    AF_allpopmax_hg38 <- rp_scores_hg38$AF_allpopmax,stringsAsFactors=FALSE),c("primer_id", "seq", "chr","hg38_pos","AF_allpopmax_hg38")))
        
         
         ##########################now hg19######################################
         
          #names(primers_hg38) <- c("primer_id","seq","chr","hg38_pos","AF_allpopmax_hg38")
         subject_hg19 <- hg19[[temp_chr]]
         products_hg19 <- matchProbePair(temp_FP,temp_RP,subject_hg19)
         amp_start_hg19 = start(products_hg19)
         fp_end_hg19 = ( start(products_hg19) + as.integer(nchar(row$FP)) - 1)
         rp_start_hg19 = ( (start(products_hg19) + width(products_hg19)) - as.integer(nchar(row$RP)) )
         amp_end_hg19 = (start(products_hg19) + width(products_hg19) - 1)
         
         
         fp_range_hg19 <- GRanges(seqnames=temp_chr, IRanges(start=amp_start_hg19:fp_end_hg19, width=1))
         fp_scores_hg19 <- gscores(mafdb_hg19,fp_range_hg19,pop=c("AF_afr","AF_amr","AF_asj","AF_eas","AF_fin","AF_nfe","AF_oth"))
         fp_scores_hg19 <- as.data.frame(fp_scores_hg19)
         
         rp_range_hg19 <- GRanges(seqnames=temp_chr, IRanges(start=rp_start_hg19:amp_end_hg19, width=1))
         rp_scores_hg19 <- gscores(mafdb_hg19,rp_range_hg19,pop=c("AF_afr","AF_amr","AF_asj","AF_eas","AF_fin","AF_nfe","AF_oth"))
         rp_scores_hg19 <- as.data.frame(rp_scores_hg19)
         
         primers_hg19 <- bind_rows(primers_hg19,setNames(data.frame(primer_id <- paste0(row$ID,"_F"),
                                                                    seq <- temp_FP,
                                                                    chr <- fp_scores_hg19$seqnames,
                                                                    hg19_pos <- fp_scores_hg19$start,
                                                                    AF_afr <- fp_scores_hg19$AF_afr,
                                                                    AF_amr <- fp_scores_hg19$AF_amr,
                                                                    AF_asj <- fp_scores_hg19$AF_asj,
                                                                    AF_eas <- fp_scores_hg19$AF_eas,
                                                                    AF_fin <- fp_scores_hg19$AF_fin,
                                                                    AF_nfe<- fp_scores_hg19$AF_nfe,
                                                                    AF_oth<- fp_scores_hg19$AF_oth,
                                                                    stringsAsFactors=FALSE),
                                                         c("primer_id", "seq", "chr","hg19_pos","AF_afr","AF_amr","AF_asj","AF_eas","AF_fin","AF_nfe","AF_oth")))
         primers_hg19 <- bind_rows(primers_hg19,setNames(data.frame(primer_id <- paste0(row$ID,"_R"),
                                                                    seq <- temp_RP,
                                                                    chr <- rp_scores_hg19$seqnames,
                                                                    hg19_pos <- rp_scores_hg19$start,
                                                                    AF_afr <- rp_scores_hg19$AF_afr,
                                                                    AF_amr <- rp_scores_hg19$AF_amr,
                                                                    AF_asj <- rp_scores_hg19$AF_asj,
                                                                    AF_eas <- rp_scores_hg19$AF_eas,
                                                                    AF_fin <- rp_scores_hg19$AF_fin,
                                                                    AF_nfe<- rp_scores_hg19$AF_nfe,
                                                                    AF_oth<- rp_scores_hg19$AF_oth,
                                                                    stringsAsFactors=FALSE),
                                                         c("primer_id", "seq", "chr","hg19_pos","AF_afr","AF_amr","AF_asj","AF_eas","AF_fin","AF_nfe","AF_oth")))
         
         
    }
 
    return(list(as.data.frame(primers_hg38), as.data.frame(primers_hg19)))
  })
     output$hg38 <- DT::renderDT(my_data()[[1]],
                                 escape=F,
                                 filter = list(position = 'bottom', clear = FALSE, plain = TRUE),
                                 rownames= F,
                                 extensions = list("ColReorder" = NULL,
                                                   "Buttons" = NULL,
                                                   "FixedColumns" = list(leftColumns=1)),
                                 options = list(
                                   dom = 'BRrltpi',
                                   autoWidth=TRUE,
                                   pageLength = 100,
                                   lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                                   ColReorder = TRUE,
                                   buttons =
                                     list(
                                       'copy',
                                       'print',
                                       list(
                                         extend = 'collection',
                                         buttons = c('csv', 'excel', 'pdf'),
                                         text = 'Download'
                                       ),
                                       I('colvis')
                                     )
                                 ))
                                 
     
    output$hg19 <- DT::renderDT(my_data()[[2]],  
                                escape=F,
                                filter = list(position = 'bottom', clear = FALSE, plain = TRUE),
                                rownames= F,
                                extensions = list("ColReorder" = NULL,
                                                  "Buttons" = NULL,
                                                  "FixedColumns" = list(leftColumns=1)),
                                options = list(
                                  dom = 'BRrltpi',
                                  autoWidth=TRUE,
                                  pageLength = 100,
                                  lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                                  ColReorder = TRUE,
                                  buttons =
                                    list(
                                      'copy',
                                      'print',
                                      list(
                                        extend = 'collection',
                                        buttons = c('csv', 'excel', 'pdf'),
                                        text = 'Download'
                                      ),
                                      I('colvis')
                                    )
                                ))
                                

  
}

shinyApp(ui, server)


