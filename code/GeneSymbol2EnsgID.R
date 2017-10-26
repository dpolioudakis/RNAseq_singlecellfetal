# Both a library and a script for mapping gene symbols (e.g. MAP1A) to their ensembl
# gene IDs (e.g. ENSG00000164024), while accounting for symbol aliases (e.g. METAP1).
# See for instance http://www.genenames.org/help/symbol-checker for motivations behind
# multi-symbol ('alias') mapping.
#
# For use as a library, the script itself can be sourced, i.e.
#    source('GeneSymbol2EnsgID.R')
#
# It can also be used as a command line utility. Use
#    Rscript GeneSymbol2EnsgID.R --help
#
# for instructions. It uses by default the most recent version of ENSEMBL for mapping
# to Ensembl gene IDs. Earlier versions can be accessed by setting the ENSEMBL_ variable
# prior to sourcing, i.e.
#    ENSEMBL_ = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", 
#                       host='sep2013.archive.ensembl.org')
#    source('GeneSymbol2EnsgID.R')

library(argparse)
library(biomaRt)
library(DBI)
library(org.Hs.eg.db)

# default statics
GOOD_CHRS_ <- c(as.character(1:22), 'X', 'Y', 'MT')

# allow user to override ENSEMBL_ version by specifying it before sourcing the utils script
if ( ! 'ENSEMBL_' %in% ls() ) {
  ENSEMBL_ <- useMart(host    = 'uswest.ensembl.org', 
                      biomart = 'ENSEMBL_MART_ENSEMBL', 
                      dataset = 'hsapiens_gene_ensembl')
}

# read arguments
get_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument('gene_list', help='The gene list', type='character')
  parser$add_argument('out_list', help='The output ensid list', type='character')
  parser$add_argument('--full_table', type='character', 
                        help='Dump the full alias table to this file')

  parser$parse_args()
}

FetchGeneAliases <- function(gene.symbol.list, verbose=F) {
  # given a list of (human) gene symbols, find aliases to them using the
  # human gene info database
  # inputs:
  #  gene.symbol.list - character vector - the desired gene symbols
  #  verbose - whether to print the number of aliases found
  #
  # returns:
  #  a dataframe containing the columns `symbol` and `alias`. If there are multiple
  #  aliases for an input gene symbol, they appear on multiple rows.

  db.con <- org.Hs.eg_dbconn()
  alias.query <- 'select * from alias, gene_info where alias._id == gene_info._id;'
  alias.symbol <- dbGetQuery(db.con, alias.query)

  gene.alias.df <- lapply(gene.symbol.list, function(sym) {
    if ( sym %in% alias.symbol$alias_symbol ) {
      gids <- subset(alias.symbol, alias_symbol == sym)$symbol
      if ( ! sym %in% gids ) {
        gids <- c(sym, gids)
      }
      data.frame('symbol'=rep(sym, length(gids)), 'alias'=gids)
    } else {
      data.frame('symbol'=sym, 'alias'=sym)
    }
  })

  gene.alias.df <- do.call(rbind, gene.alias.df)
  
  if (verbose) {
    print(paste('Aliases found:', length(setdiff(gene.alias.df$alias, gene.alias.df$symbol))))
  }

  gene.alias.df
}

LookupEnsemblIDs <- function(gene.frame, chr.filt=GOOD_CHRS_, verbose=F, mart=ENSEMBL_) {
  # look up gene ensembl ids from gene symbols, using ensembl.org
  # inputs:
  #   gene.frame - data frame containing the columns `symbol` and `alias` (these can be the same)
  #   chr.filt - a white-list of chromosomes from which to accept genes. Default: 1-22, X, and Y
  #   verbose - whether to print information on how many genes could be mapped
  #   mart - the biomaRt to use, if different from the default ensembl
  #
  # returns:
  #   a data frame with the columns `symbol`, `alias`, and `ensid`.
  #   Multiple ENSIDs for the same alias will are placed in separate rows. 

  bmr <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name'), 
               filters    = c('external_gene_name', 'chromosome_name'), 
               values     = list(gene.frame$alias, chr.filt), mart=ENSEMBL_)

  gene.ensembl.df <- lapply(gene.frame$alias, function(galias) {
    if ( galias %in% bmr$external_gene_name ) {
      ensids <- subset(bmr, external_gene_name == galias)$ensembl_gene_id
      symbol <- subset(gene.frame, alias==galias)$symbol
      data.frame('symbol'=rep(symbol, length(ensids)), 'alias'=rep(galias, length(ensids)), 
                   'ensid'=ensids)
    } else {
      symbol <- subset(gene.frame, alias==galias)$symbol
      data.frame('symbol'=symbol, 'alias'=galias, 'ensid'='NA')
    }
  })

  gene.ensembl.df <- unique(do.call(rbind, gene.ensembl.df))

  if ( verbose ) {
    # for convenience rename
    gene.symbol.list <- gene.frame$symbol
    print(paste('Starting symbol size:', length(unique(gene.symbol.list))))
    gene.ensid.size <- length(unique(gene.ensembl.df$ensid[gene.ensembl.df$ensid != 'NA']))
    print(paste('Resulting ensid list size:', gene.ensid.size))

    # check for genes which failed to convert and display them
    failed <- sapply(gene.symbol.list, function(sym) {
      ssub <- subset(gene.ensembl.df, symbol==sym)
      all(ssub$ensid == 'NA')
    })

    if ( any(failed) ) {
      print('Failed to convert:')
      print(subset(gene.ensembl.df, symbol %in% gene.symbol.list[failed]))
    }

  }

  gene.ensembl.df
}

main <- function() {
  args <- get_args()
  gene.symbol.list <- read.table(args$gene_list, header=F, stringsAsFactors=F)[,1]
  alias.df <- FetchGeneAliases(gene.symbol.list, verbose=T)
  good.chrs <- c(GOOD_CHRS_, 'CHR_HSCHR6_MHC_SSTO_CTG1')  # HLA-DRB4
  gene.ensembl.df <- LookupEnsemblIDs(alias.df, chr.filt=good.chrs, verbose=T)
    
  # write the converted ids
  cat(gene.ensembl.df$ensid[gene.ensembl.df$ensid != 'NA'], file=args$out_list, sep='\n')

  # if the full table was requested, write it out
  if ( ! is.null(args$full_table) ) {
    write.table(gene.ensembl.df, quote=F, row.names=F, file=args$full_table)
  }
}

if ( ! interactive() ) {
  main()
}
