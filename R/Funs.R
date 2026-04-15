#' Fully Parameterizable ANOVA test
#'
#' @description Function allowing the user to perform ANOVA tests (in a strict sense: one continuous and normally-distributed response variable and one or more factorial/categorical explaining variable(s)),
#'  with the possibility to specify the type of sum of squares (1, 2 or 3), the types of factors (Fixed or Random) and their relationships (crossed or nested).
#'  This function is built on the '?rstatix::anova_test()'.The resulting ouptuts are the same as in 'SAS' software.
#'
#' @param InputData The name of the input dataset.
#' @param form To be written as _Response variable ~ Explaining variable(s)_. If there are multiple _Explaining variables_, they should be separated by a _*_.
#' @param variables_type Type of factors in a vector. One argument per _Explaining variable_, in the same order as in the formula: 'Fxd' for Fixed factors, or 'Rndm' for Random factors
#' @param variables_relationships Relationships between factors in a vector. One argument per _Explaining variable_, in the same order as in the formula:
#'  _NA_ if the factor is not nested in any other factor, or the name of the factor it is nested into.
#' @param sig_threshold The threshold below which a p-value is considered significant. By default, it is set to the usual value of 0.05.
#' @param sum_of_squares Type of sum of squares used in the test: 1, 2 or 3. These numbers correspond to those in 'SAS' software. By default, the type 3 is used as it is the most often recommended one.
#' @param print_ges If _TRUE_, prints effect sizes (ges = generalized eta squared) in the result table. See '?rstatix::anova_test()' for more details. It is _FALSE_ by default.
#' @param output_table_only If _TRUE_, reports only the result table and not the text below specifying the parameters used for the test, allowing for the output table to be directly usable or stored in an object for example. It is _FALSE_ by default.
#'
#' @return
#' ANOVA result table.
#'
#' There is one row per effect tested and a last row for the residuals.
#'
#' The first column ("Effect") specifies the effects tested.
#'
#' The second column ("SS") gives the Sum of Squares associated with the effects tested.
#'
#' The third column ("DF") gives the Degrees of Freedom associated with the effects tested.
#'
#' The fourth column ("MSS") gives the Mean Sum of Squares (=SS/DF) associated with the effects tested.
#'
#' The fifth column ("Fval") gives the F statistic value associated with the effects tested.
#'
#' The sixth column ("pval") gives the p-value associated with the effects tested.
#'
#' The seventh column ("Sig") shows a star if the associated p-value is below the statistical threshold of 0.05 by default or the one indicated during the call of the function for each effects tested.
#'
#' The optional eighth column ("ges") gives the effect sizes (ges = generalized eta squared) associated with the effects tested.
#'
#' If the argument _output_table_only_ is _FALSE_ (which is the case by default), the function also specifies below this result table which types of Sum of Squares were used,
#' what is the response variable, which factors were fixed or random, and if there are at least two factors, it specifies the relationships between these factors.
#'
#' You will find warnings below if the setup is not appropriate.
#'
#' @export
#'
#' @examples
#'data("Butterfly")
#'# Example of ANOVA I
#'FullyParamANOVA(InputData=Butterfly, form=Wingspan~Site, variables_type = c('Fxd'))
#'
#'# Example of ANOVA II
#'FullyParamANOVA(InputData=Butterfly, form=Wingspan~Patch*Site,
#'variables_type = c('Rndm','Fxd'), variables_relationships = c("Site",NA))
#'
#'# Example of ANOVA III
#'FullyParamANOVA(InputData=Butterfly, form=Wingspan~Sex*Patch*Site,
#'variables_type = c('Fxd','Rndm','Fxd'), variables_relationships = c(NA,"Site",NA))
#'
#'# Example of ANOVA III, using arguments different than those by default
#'# (changing the threshold below which a p-value is considered significant to 0.1,
#'# the types of sum of squares to 2, adding the effect sizes (ges = generalized eta squared)
#'# in the output table, and outputting only the table to be able to use the information within
#'# or store it in an object)
#'FullyParamANOVA(InputData=Butterfly, form=Wingspan~Sex*Patch*Site, sig_threshold=0.1,
#'sum_of_squares = 2, variables_type = c('Fxd','Rndm','Fxd'),
#'variables_relationships = c(NA,"Site",NA), print_ges = TRUE, output_table_only = TRUE)
#'
#' @importFrom dplyr mutate filter n summarise
#' @importFrom mlr3misc extract_vars
#' @importFrom rstatix anova_test
#' @importFrom tibble is_tibble
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
FullyParamANOVA<- function(InputData,form,variables_type=NULL,variables_relationships= c(NA),sig_threshold=0.05,sum_of_squares=3,print_ges=FALSE,output_table_only=FALSE){
  # Check if inputs are correct
  # The threshold value below which a p-value is considered significant must be a number between 0 and 1
  if (!is.numeric(sig_threshold) | sig_threshold<=0 | sig_threshold>=1){
    stop("The threshold value below which a p-value is considered significant must be a number between 0 and 1",call. = FALSE)
  }
  # The type of sum of squares must be 1, 2, or 3
  if (!is.numeric(sum_of_squares) | !sum_of_squares %in% c(1,2,3)){
    stop("The type of sum of squares must be 1, 2, or 3",call. = FALSE)
  }
  # The argument to indicate if ges must be printed must be logical (i.e. TRUE or FALSE)
  if (!is.logical(print_ges)){
    stop("The argument to indicate if ges must be printed must be logical (i.e. TRUE or FALSE)",call. = FALSE)
  }
  # The argument to indicate if only the output table must be printed must be logical (i.e. TRUE or FALSE)
  if (!is.logical(output_table_only)){
    stop("The argument to indicate if only the output table must be printed must be logical (i.e. TRUE or FALSE)",call. = FALSE)
  }
  # Dataset must be a data frame or a tibble
  if (!is.data.frame(InputData) & !tibble::is_tibble(InputData)) {
    stop("Please use a correct dataframe or tibble as input.",call. = FALSE)
  }
  # Works better with a data frame than a tibble
  if (tibble::is_tibble(InputData)){
    InputData<-as.data.frame(InputData)
  }
  # Variable must not have ":" in their names
  if (sum(grepl(":",mlr3misc::extract_vars(form)$rhs,fixed = TRUE))>0){
    stop("Variable must not have ":" in their names",call. = FALSE)
  }
  # Explanatory variables must be separated by a '*' in the formula
  if (sum(grepl("+",as.character(form),fixed = TRUE)>0)){
    stop("Explanatory variables must only be separated by a '*' in the formula",call. = FALSE)
  }
  # There must not be more than 3 factors
  if (length(mlr3misc::extract_vars(form)$rhs)>3) {
    stop("This function does not support models with more than 3 explanatory variables.",call. = FALSE)
  }
  #  There must be only one continuous response variable
  if (length(mlr3misc::extract_vars(form)$lhs)>1 | !is.numeric(InputData[,mlr3misc::extract_vars(form)$lhs])){
    stop("There must be only one continuous numerical response variable",call. = FALSE)
  }
  # The types of factors must be mentioned, there must be one type element per factor, and these elements must either be "Fxd" for fixed or "Rndm" for random
  if (is.null(variables_type) |
      (length(mlr3misc::extract_vars(form)$rhs))!=length(variables_type) |
      sum(variables_type %in% c("Fxd","Rndm"))!=length(variables_type)) {
    stop("Please mention the factors types ('Fxd' for Fixed and 'Rndm' for Random) for each factor, in the same order as in the formula.\nExample for a model with two factors, the first one fixed, the second random: variables_type=c('Fxd','Rndm')",call. = FALSE)
  }
  # The relationships between factors must be mentioned, there must be one relationship element per factor, and these elements must either be "NA" if the factor is not nested, or the name of the factor or the interaction it is nested into
  if (is.null(variables_relationships)|
      (length(mlr3misc::extract_vars(form)$rhs))!=length(variables_relationships) |
      sum(variables_relationships %in% c(NA,
                                         mlr3misc::extract_vars(form)$rhs,
                                         (expand.grid(Var1=mlr3misc::extract_vars(form)$rhs,Var2=mlr3misc::extract_vars(form)$rhs) %>% dplyr::filter(.data$Var1!=.data$Var2) %>% dplyr::mutate(Interaction=base::paste(.data$Var1,.data$Var2,sep = ":")))$Interaction))!=length(variables_relationships) |
      sum(variables_relationships[which(!is.na(variables_relationships))] %in% c(expand.grid(mlr3misc::extract_vars(form)$rhs,mlr3misc::extract_vars(form)$rhs) %>% dplyr::filter(.data$Var1!=.data$Var2) %>% dplyr::mutate(Interaction=base::paste(.data$Var1,.data$Var2,sep = ":")))$Interaction)>1) {
    stop("Please write the relationships between factors (one per factor, in the same order).\nExample 1 for a model with two factors, the second one nested in the first one (called X1 here): variables_relationships=c(NA,'X1')\nExample 2 if the third factor is nested in the combination of the two first others, which are crossed (X1 and X2 here): variables_relationships=c(NA,NA,'X1:X2')\nExample 3 if the third factor is nested in the second factor (X2 here) which is also nested in the first factor (X1 here): variables_relationships=c(NA,'X1','X1:X2')",call. = FALSE)
  }
  # The relationships must be logical (at least one non nested factor and no self nesting (even in interaction), no nesting cycle and correct order when full nested design)
  if (sum(is.na(variables_relationships))==0) {
    stop("Relationships between factors are not logical.\nThere should be at least one factor that is not nested to anything.\nPlease check again.",call. = FALSE)
  }
  if (sum(!is.na(variables_relationships))>0){
    for (i in 1:length(variables_relationships[which(!is.na(variables_relationships))])){
      if(variables_relationships[which(!is.na(variables_relationships))][i]==mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][i]){
        stop("Relationships between factors are not logical.\nA factor can not be nested in itself.\nPlease check again.",call. = FALSE)
      }
    }
    if (sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))>0){
      if (grepl(mlr3misc::extract_vars(form)$rhs[which(grepl(":",variables_relationships,fixed = TRUE))],variables_relationships[which(grepl(":",variables_relationships,fixed = TRUE))],fixed = TRUE)){
        stop("Relationships between factors are not logical.\nA factor can not be nested in itself (here in the interaction).\nPlease check again.",call. = FALSE)
      }
    }
    if (sum(!is.na(variables_relationships))==2 & sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==1){
      if(variables_relationships[which(!grepl(":",variables_relationships,fixed=TRUE) & !is.na(variables_relationships))]!=mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))] |
         !grepl(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],variables_relationships[which(grepl(":",variables_relationships,fixed=TRUE))],fixed=TRUE) |
         !grepl(mlr3misc::extract_vars(form)$rhs[which(!grepl(":",variables_relationships,fixed=TRUE) & !is.na(variables_relationships))],variables_relationships[which(grepl(":",variables_relationships,fixed=TRUE))],fixed=TRUE)){
        stop("Relationships between factors are not logical in this fully nested design.\nPlease check again.",call. = FALSE)
      }
    }
    if (sum(!is.na(variables_relationships))==2){
      if (sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==0 &
          length(unique(variables_relationships[which(!is.na(variables_relationships))]))!=1){
        if(variables_relationships[which(!is.na(variables_relationships))][1]==mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2] &
           variables_relationships[which(!is.na(variables_relationships))][2]==mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1]){
          stop("Relationships between factors are not logical.\nThere should not be any factor nesting loop (eg A nested in B and B nested in A).\nPlease check again.",call. = FALSE)
        }
        else stop("Relationships between factors are not logical or not correctly written.\nIf you tried to write a fully nested design, you should follow this example: if the third factor is nested in the second factor (X2 here) which is also nested in the first factor (X1 here): variables_relationships=c(NA,'X1','X1:X2')\nIf you did not try to do such a design, then please check again.",call. = FALSE)
      }
    }
  }
  # If only one factor and only one observation per level, then nothing can be tested
  if (((length(mlr3misc::extract_vars(form)$rhs)==1 & length(variables_type)==1)) &
      (sum((InputData %>% dplyr::group_by(dplyr::across(dplyr::all_of(mlr3misc::extract_vars(form)$rhs))) %>% dplyr::summarise(n=dplyr::n(),.groups = "keep"))$n>1)==0)) {
    stop("If there is only 1 observation then nothing can be tested",call. = FALSE)
  }
  # The nested factors should be random
  if (sum(variables_type=="Rndm" & !is.na(variables_relationships))!=sum(!is.na(variables_relationships))){
    stop("Nested factors are by definition random factors.\nPlease adapt factors type argument.",call. = FALSE)
  }
  #  The factors must be factors, if not, transformed into factors and warned
  if (length(mlr3misc::extract_vars(form)$rhs)==1 & !is.factor(InputData[,mlr3misc::extract_vars(form)$rhs])){
    warning(base::paste0("The function transformed ",
                         mlr3misc::extract_vars(form)$rhs,
                         " as factor because it was not the case yet.\nTo remove this warning, please transform ",
                         mlr3misc::extract_vars(form)$rhs,
                         " into a factor in your input dataset, or check if you put the correct factors in the formula."),call. = FALSE)
    InputData[,mlr3misc::extract_vars(form)$rhs]<-as.factor(InputData[,mlr3misc::extract_vars(form)$rhs])
  }
  else if (length(mlr3misc::extract_vars(form)$rhs)>1 & sum(sapply(InputData[,mlr3misc::extract_vars(form)$rhs],is.factor))!=length(mlr3misc::extract_vars(form)$rhs)){
    warning(base::paste0("The function transformed ",
                         base::paste0(names(sapply(InputData[,mlr3misc::extract_vars(form)$rhs],is.factor))[which(!sapply(InputData[,mlr3misc::extract_vars(form)$rhs],is.factor))],collapse = " & "),
                         " as factor because it was not the case yet.\nTo remove this warning, please transform ",
                         base::paste0(names(sapply(InputData[,mlr3misc::extract_vars(form)$rhs],is.factor))[which(!sapply(InputData[,mlr3misc::extract_vars(form)$rhs],is.factor))],collapse = " & "),
                         " into a factor in your input dataset, or check if you put the correct factors in the formula."),call. = FALSE)
    InputData[,mlr3misc::extract_vars(form)$rhs] <- lapply(InputData[,mlr3misc::extract_vars(form)$rhs],as.factor)
  }
  # Warning if the data is unbalanced
  if (length(unique((InputData %>% dplyr::group_by(dplyr::across(dplyr::all_of(mlr3misc::extract_vars(form)$rhs))) %>% dplyr::summarise(n=dplyr::n(),.groups = "keep"))$n))>1){
    warning("Data is unbalanced. Please take these results cautiously.",call. = FALSE)
  }




  # Compute results
  # If there is at least two observations in at least one combination
  if (sum((InputData %>% dplyr::group_by(dplyr::across(dplyr::all_of(mlr3misc::extract_vars(form)$rhs))) %>% dplyr::summarise(n=dplyr::n(),.groups = "keep"))$n>1)>0){

    # If all fixed and crossed, or has just one factor, output from anova_test() is correct (if there is more than 1 observation per combination)
    if ((length(variables_type)==sum(variables_type=='Fxd') & length(variables_relationships)==sum(is.na(variables_relationships))) |
        (length(mlr3misc::extract_vars(form)$rhs)==1 & length(variables_type)==1)) {
      Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")
    }

    # If not, has to change the computations
    # Cases when there are two factors
    else if (length(mlr3misc::extract_vars(form)$rhs)==2 & length(variables_type)==2){
      # If there are no nested factors, but one is random
      if (length(variables_relationships)==sum(is.na(variables_relationships)) & sum(variables_type=='Rndm')==1){
        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)
      }

      # If there are no nested factors and all random
      else if (length(variables_relationships)==sum(is.na(variables_relationships)) & sum(variables_type=='Rndm')==2){
        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)
      }

      # If one factor is nested in the other (fixed or random does not change the results)
      else if (length(variables_relationships)>sum(is.na(variables_relationships))){
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":")))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[!is.na(variables_relationships)],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[!is.na(variables_relationships)],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[!is.na(variables_relationships)],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[!is.na(variables_relationships)],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[!is.na(variables_relationships)],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[!is.na(variables_relationships)],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],"(",mlr3misc::extract_vars(form)$rhs[is.na(variables_relationships)],")")
      }
    }


    # Cases when there are three factors
    else if (length(mlr3misc::extract_vars(form)$rhs)==3 & length(variables_type)==3){
      # If there are no nested factors, but one is random
      if (length(variables_relationships)==sum(is.na(variables_relationships)) & sum(variables_type=='Rndm')==1){
        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                     c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                       mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1]))],collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                       c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                         mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1]))],collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                              c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                                mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1]))],collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                     c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                       mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2]))],collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                       c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                         mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2]))],collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                              c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                                mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2]))],collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],collapse = ":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],collapse = ":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],collapse = ":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],collapse = ":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],collapse = ":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)
      }

      # If there are no nested factors, but two are random
      else if (length(variables_relationships)==sum(is.na(variables_relationships)) & sum(variables_type=='Rndm')==2){
        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$Fval<-NA
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$p<-NA

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                  c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],
                                                                                                                    mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                     c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],
                                                                                                                       mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                       c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],
                                                                                                                         mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                  c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],
                                                                                                                    mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                              c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],
                                                                                                                                mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                              c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],
                                                                                                                                mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                  c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],
                                                                                                                    mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                     c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],
                                                                                                                       mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                       c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],
                                                                                                                         mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                  c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],
                                                                                                                    mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                              c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],
                                                                                                                                mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                              c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],
                                                                                                                                mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")]))],collapse = ":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)
      }

      # If there are no nested factors, but they are all random
      else if (length(variables_relationships)==sum(is.na(variables_relationships)) & sum(variables_type=='Rndm')==3){
        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$Fval<-NA
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$p<-NA

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$Fval<-NA
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$p<-NA

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[3],]$Fval<-NA
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[3],]$p<-NA

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[2]),collapse = ":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[2]),collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[2]),collapse = ":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[2]),collapse = ":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[2]),collapse = ":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[2]),collapse = ":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[1],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[2],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[2],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[2],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[2],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[2],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(c(mlr3misc::extract_vars(form)$rhs[2],mlr3misc::extract_vars(form)$rhs[3]),collapse = ":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs,collapse = ":"),]$DFn, lower.tail = FALSE)
      }

      # If nested model of type I (fully nested; fixed or random does not change the results)
      else if (sum(!is.na(variables_relationships))==2 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==1){
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & grepl(":",variables_relationships,fixed = TRUE))],sep=":")))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$DFn, lower.tail = FALSE)


        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],"(",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],")")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & grepl(":",variables_relationships,fixed = TRUE))],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & grepl(":",variables_relationships,fixed = TRUE))],"(",
                       mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))],"(",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],"))")
      }

      # If nested model of type II (fixed or random does not change the results)
      else if(sum(!is.na(variables_relationships))==2 & sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==0){
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":")))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$Fval<-NA
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$p<-NA

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$DFn, lower.tail = FALSE)


        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$Fval<-
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$DFn) /
          (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                                           mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$SSn /
             Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$DFn)
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$p<-
          stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$Fval,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$DFn,
                    Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                                                    mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$DFn, lower.tail = FALSE)

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],"(",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],")")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],"(",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],")")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$Effect<-
          base::paste(base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],"(",
                                   mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],")"),
                      base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],"(",
                                   mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],")"),sep=":")
      }

      # If nested model of type III
      else if (sum(!is.na(variables_relationships))==1 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==1){
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],"+",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":")))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")
        # If both factors that are not nested are fixed
        if (sum(variables_type=="Fxd" & is.na(variables_relationships))==2){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                          mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                          mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
        }

        # If one of the factors that are not nested is fixed and the other random
        else if (sum(variables_type=="Fxd" & is.na(variables_relationships))==1){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Fxd")],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Fxd")],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Fxd")],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Fxd")],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Fxd")],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Fxd")],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Rndm")],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Rndm")],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Rndm")],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Rndm")],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Rndm")],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & variables_type=="Rndm")],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                          mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                          mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
        }

        # If both factors that are not nested are random
        else if (sum(variables_type=="Fxd" & is.na(variables_relationships))==0){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                          mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                             mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                               mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                          mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

        }

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],
                                                                        mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],"(",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1],":",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2],")")
      }

      # If nested model of type IV
      else if (sum(!is.na(variables_relationships))==1 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==0) {
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],"+",
                                                 variables_relationships[which(!is.na(variables_relationships))],"+",
                                                 base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":")))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")
        # If both the crossed factor and the nesting factor are fixed
        if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Fxd" &
            variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Fxd"){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                            mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
        }

        # If the crossed factor is fixed and the nesting factor is random
        else if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Fxd" &
                 variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Rndm"){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                            mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
        }

        # If the crossed factor is random and the nesting factor is fixed
        else if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Rndm" &
                 variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Fxd"){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                            mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-NA

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                          mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                          mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
        }

        # If the crossed factor and the nesting factor are random
        else if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Rndm" &
                 variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Rndm"){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                            mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-NA

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                          variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                          mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                          mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
        }
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],"(",
                       variables_relationships[which(!is.na(variables_relationships))],")")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                 mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                        variables_relationships[which(!is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],":",
                       mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],"(",
                       variables_relationships[which(!is.na(variables_relationships))],")")
      }
    }
  }

  # If there is one observation per combination of levels of factors
  # If expect no highest interaction effect, can consider highest interaction as a measure of residuals for tests against CMR
  else {
    # Cases when there are two factors
    if (length(mlr3misc::extract_vars(form)$rhs)==2 & length(variables_type)==2){
      # If both are crossed
      if(length(variables_relationships)==sum(is.na(variables_relationships))){
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[1],"+",
                                                 mlr3misc::extract_vars(form)$rhs[2]))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")
        # If all fixed
        if(length(variables_type)==sum(variables_type=='Fxd')){
          warning("Only one observation per combination of levels of factors, so the interaction effect cannot be tested. If the assumption of the non existence of the interaction is true, then the effects of both factors can be tested using the interaction as a measure of residual variation, which gives the results above. If the assumption of the non existence of the interaction is false, then these results would be biased.",call. = FALSE)
        }
        # If one fixed and one random
        else if(sum(variables_type=='Fxd')==1 & sum(variables_type=='Rndm')==1){
          warning("Only one observation per combination of levels of factors, so the interaction effect cannot be tested. If the assumption of the non existence of the interaction is true, then the effect of the random factor can be tested using the interaction as a measure of residual variation, which gives the results above. If the assumption of the non existence of the interaction is false, then the test for the random factor would be biased. In either case, the test for the fixed factor is not biased",call. = FALSE)
        }
        # If all random
        else if(length(variables_type)==sum(variables_type=='Rndm')){
          warning("Only one observation per combination of levels of factors, so the interaction effect cannot be tested. The tests for both factors are however not biased.",call. = FALSE)
        }
      }
      # If one nested within the other
      else if(sum(is.na(variables_relationships))==1){
        stop("You have considered a nested factor that should be considered as a replicate (nested factor with only one observation per level). Please remove this nested factor from your formula, and analyze your data with an ANOVA I design.",call. = FALSE)
      }
    }
    # Cases when there are three factors
    if (length(mlr3misc::extract_vars(form)$rhs)==3 & length(variables_type)==3){
      # If all crossed
      if(length(variables_relationships)==sum(is.na(variables_relationships))){
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[1],"+",
                                                 mlr3misc::extract_vars(form)$rhs[2],"+",
                                                 mlr3misc::extract_vars(form)$rhs[3],"+",
                                                 mlr3misc::extract_vars(form)$rhs[1],":",mlr3misc::extract_vars(form)$rhs[2],
                                                 mlr3misc::extract_vars(form)$rhs[1],":",mlr3misc::extract_vars(form)$rhs[3],
                                                 mlr3misc::extract_vars(form)$rhs[2],":",mlr3misc::extract_vars(form)$rhs[3]))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")
        #If all fixed
        if(length(variables_type)==sum(variables_type=='Fxd')){
          warning("Only one observation per combination of levels of factors, so the triple interaction effect cannot be tested. If the assumption of the non existence of this triple interaction is true, then the effects of all three factors and their double interaction can be tested using the triple interaction as a measure of residual variation, which gives the results above. If the assumption of the non existence of the triple interaction is false, then these results would be biased.",call. = FALSE)
        }
        # If two fixed and one random
        else if(sum(variables_type=='Fxd')==2 & sum(variables_type=='Rndm')==1){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                       c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                         mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1]))],collapse = ":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                         c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                           mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1]))],collapse = ":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                                c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                                  mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][1]))],collapse = ":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                       c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                         mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2]))],collapse = ":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                         c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                           mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2]))],collapse = ":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(mlr3misc::extract_vars(form)$rhs %in%
                                                                                                                                c(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],
                                                                                                                                  mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")][2]))],collapse = ":"),]$DFn, lower.tail = FALSE)
          warning("Only one observation per combination of levels of factors, so the triple interaction effect cannot be tested. If the assumption of the non existence of the triple interaction is true, then the three effects related to the random factor (i.e. the effect of the random factor, and the two double interactions involved with the random factor) can be tested using the triple interaction as a measure of residual variation, which gives the results above. If the assumption of the non existence of the triple interaction is false, then the tests for the three effects related to the random factor would be biased. In either case, the tests for both fixed factors, as well as the interaction between those two are not biased",call. = FALSE)
        }
        # If one fixed and two random
        else if(sum(variables_type=='Fxd')==1 & sum(variables_type=='Rndm')==2){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],]$p<-NA

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][1],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")][2],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste0(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ":"),]$DFn, lower.tail = FALSE)
          warning("Only one observation per combination of levels of factors, so the triple interaction effect cannot be tested. If the assumption of the non existence of the triple interaction is true, then the effect of the double interaction between the two random factors can be tested using the triple interaction as a measure of residual variation, which gives the results above. If the assumption of the non existence of the triple interaction is false, then this test of the interaction between the two random factors would be biased. In either case, the four other computable tests (each of the two random factors and the two double interactions between the fixed factor and either one of the random factors) are not biased.",call. = FALSE)
        }
        # If all random
        else if(sum(variables_type=='Rndm')==3){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[1],]$p<-NA

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[2],]$p<-NA

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[3],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[3],]$p<-NA
          warning("Only one observation per combination of levels of factors, so the triple interaction effect cannot be tested. The three computable tests (all double interactions) are however not biased.",call. = FALSE)
        }
      }
      # If nested model of type I (fully nested; fixed or random does not change the results)
      else if (sum(!is.na(variables_relationships))==2 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==1){
        stop("You have considered a nested factor that should be considered as a replicate (nested factor with only one observation per level, i.e. the nested factor nested within another nested factor). Please remove this nested factor from your formula, and analyze your data with an nested ANOVA II design.",call. = FALSE)
      }
      # If nested model of type II (fixed or random does not change the results)
      else if(sum(!is.na(variables_relationships))==2 & sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==0){
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":")))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$Fval<-NA
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],]$p<-NA

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1],"(",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],")")

        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2],"(",
                       mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))],")")
        warning("Only one observation per combination of levels of factors, so the effect of the interaction between the two nested factors cannot be tested. The two computable tests (both nested factors) are however not biased.",call. = FALSE)
      }

      # If nested model of type III
      else if (sum(!is.na(variables_relationships))==1 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==1){
        stop("You have considered a nested factor that should be considered as a replicate (nested factor with only one observation per level, i.e. the nested factor nested within the combination of the other two factors). Please remove this nested factor from your formula, and analyze your data with an crossed ANOVA II design.",call. = FALSE)
      }
      # If nested model of type IV
      else if (sum(!is.na(variables_relationships))==1 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==0) {
        new_form<-stats::as.formula(base::paste0(mlr3misc::extract_vars(form)$lhs,"~",
                                                 mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],"+",
                                                 variables_relationships[which(!is.na(variables_relationships))],"+",
                                                 base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),"+",
                                                 base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":")))

        Results_of_ANOVA_test<-data.frame(rstatix::anova_test(InputData,new_form,type=sum_of_squares,detailed = TRUE)) %>% dplyr::rename("Fval"="F")
        # If both the crossed factor and the nesting factor are fixed
        if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Fxd" &
            variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Fxd"){


          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
          warning("Only one observation per combination of levels of factors, so the effect of the interaction between the crossed factor and the nested factor cannot be tested. If the assumption of the non existence of the interaction between the crossed factor and the nested factor is true, then the effect of the nested factor can be tested using the interaction between the crossed factor and the nested factor as a measure of residual variation, which gives the results above. If the assumption of the non existence of the interaction between the crossed factor and the nested factor is false, then this test of the nested factor would be biased. In either case, the three other computable tests (those that do not involve the nested factor) are not biased",call. = FALSE)
        }

        # If the crossed factor is fixed and the nesting factor is random
        else if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Fxd" &
                 variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Rndm"){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                            mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                             mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                               mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                                      mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)
          warning("Only one observation per combination of levels of factors, so the effect of the interaction between the crossed factor and the nested factor cannot be tested. If the assumption of the non existence of the interaction between the crossed factor and the nested factor is true, then the effect of the nested factor can be tested using the interaction between the crossed factor and the nested factor as a measure of residual variation, which gives the results above. If the assumption of the non existence of the interaction between the crossed factor and the nested factor is false, then this test of the nested factor would be biased. In either case, the three other computable tests (those that do not involve the nested factor) are not biased",call. = FALSE)
        }

        # If the crossed factor is random and the nesting factor is fixed
        else if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Rndm" &
                 variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Fxd"){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-NA
          warning("Only one observation per combination of levels of factors, so the effect of the interaction between the crossed factor and the nested factor cannot be tested. The three other computable tests (the crossed factor, the interaction between the crossed factor and the nesting factor, as well as the nested factor) are however not biased",call. = FALSE)
        }

        # If the crossed factor and the nesting factor are random
        else if (variables_type[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]=="Rndm" &
                 variables_type[which(mlr3misc::extract_vars(form)$rhs==variables_relationships[which(!is.na(variables_relationships))])]=="Rndm"){
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval<-
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                          mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                            mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn) /
            (Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                      mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                             variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$SSn /
               Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                        mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                               variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn)
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                       mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$p<-
            stats::pf(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$Fval,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                   mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],]$DFn,
                      Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) &
                                                                                                                               mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])],
                                                                                      variables_relationships[which(!is.na(variables_relationships))],sep=":"),]$DFn, lower.tail = FALSE)

          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$Fval<-NA
          Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==variables_relationships[which(!is.na(variables_relationships))],]$p<-NA
          warning("Only one observation per combination of levels of factors, so the effect of the interaction between the crossed factor and the nested factor cannot be tested. The three other computable tests (the crossed factor, the interaction between the crossed factor and the nesting factor, as well as the nested factor) are however not biased",call. = FALSE)
        }
        Results_of_ANOVA_test[Results_of_ANOVA_test$Effect==base::paste(variables_relationships[which(!is.na(variables_relationships))],
                                                                        mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],sep=":"),]$Effect<-
          base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))],"(",
                       variables_relationships[which(!is.na(variables_relationships))],")")
      }
    }
  }




  if (!exists("Results_of_ANOVA_test")) stop("Case not tackled yet. Please contact the author of the function.",call. = FALSE)

  # If there are some tests that were not computable, tell the user
  if (sum(is.na(Results_of_ANOVA_test))>0){
    warning("Some tests were not applicable with this experimental design, leading to NAs in F statistics values and pvalues of these tests.\nPlease refer to appropriate statistical resources for a better understanding of this particular result.",call. = FALSE)
  }

  # Print results
  if (output_table_only==FALSE){
    # Print results table
    ifelse(print_ges==TRUE,ResultsToPrint<-(data.frame(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect!="(Intercept)",]) %>%
                                              dplyr::rename(Sig=.data$p..05,SS=.data$SSn,DF=.data$DFn,pval=.data$p) %>%
                                              dplyr::mutate(Sig=ifelse(.data$pval<sig_threshold,"*",""),
                                                            MSS=formatC(signif(.data$SS/.data$DF, digits=5)),
                                                            Fval=formatC(signif(.data$Fval, digits=5)),
                                                            pval=formatC(signif(.data$pval, digits = 5))) %>%
                                              base::rbind(c("Residuals",Results_of_ANOVA_test$SSd[1],"",Results_of_ANOVA_test$DFd[1],"","","","","",formatC(signif(unique(Results_of_ANOVA_test$SSd/Results_of_ANOVA_test$DFd),digits = 5)))) %>%
                                              dplyr::mutate(SSd=NULL,DFd=NULL))[,c("Effect","SS","DF","MSS","Fval","pval","Sig","ges")],
           ResultsToPrint<-(data.frame(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect!="(Intercept)",]) %>%
                              dplyr::rename(Sig=.data$p..05,SS=.data$SSn,DF=.data$DFn,pval=.data$p) %>%
                              dplyr::mutate(Sig=ifelse(.data$pval<sig_threshold,"*",""),
                                            MSS=formatC(signif(.data$SS/.data$DF, digits=5)),
                                            Fval=formatC(signif(.data$Fval, digits=5)),
                                            pval=formatC(signif(.data$pval, digits = 5))) %>%
                              base::rbind(c("Residuals",Results_of_ANOVA_test$SSd[1],"",Results_of_ANOVA_test$DFd[1],"","","","","",formatC(signif(unique(Results_of_ANOVA_test$SSd/Results_of_ANOVA_test$DFd),digits = 5)))) %>%
                              dplyr::mutate(SSd=NULL,DFd=NULL))[,c("Effect","SS","DF","MSS","Fval","pval","Sig")])
    print(ResultsToPrint,row.names=FALSE)

    # Add wich type of sum of square was used
    cat(base::paste0("\nANOVA test with Type ",sum_of_squares, " Sum of Squares, and a significance threshold of ", sig_threshold))

    # Mention what was the response variable
    cat(base::paste0("\nResponse variable: ",mlr3misc::extract_vars(Wingspan~Sex)$lhs))

    # Mention what were the types of the factors involved
    cat("\nFactors types: ")
    if (length(mlr3misc::extract_vars(form)$rhs)==1 & length(variables_type)==1){
      cat(base::paste0(mlr3misc::extract_vars(form)$rhs," is ",ifelse(variables_type[1]=="Fxd","Fixed","Random")))
    }
    else if (sum(variables_type=="Fxd")==length(variables_type)){
      cat("All factors are Fixed")
    }
    else if (sum(variables_type=="Rndm")==length(variables_type)){
      cat("All factors are Random")
    }
    else {
      cat(base::paste0("Fixed factors = ",base::paste(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Fxd")],collapse = ",")))
      cat(base::paste0(" | Random factors = ",base::paste(mlr3misc::extract_vars(form)$rhs[which(variables_type=="Rndm")],collapse = ",")))
    }

    # If there are more than one factor, mention the relationships between the factors involved
    if (length(variables_relationships)>=2){
      cat("\nRelationships between factors: ")

      #If there are two factors
      if (length(variables_relationships)==2){
        if (sum(is.na(variables_relationships))==length(variables_relationships)){
          cat("All factors are crossed")
        }
      else {
          cat(base::paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))]," is nested within ", variables_relationships[which(!is.na(variables_relationships))]))
      }
    }

    #If there are three factors
    else if (length(variables_relationships)==3){
      if (sum(is.na(variables_relationships))==length(variables_relationships)){
        cat("All factors are crossed")
      }
      else if (sum(!is.na(variables_relationships))==2 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==1){
        cat(paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))], " is nested within ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))]," | ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & grepl(":",variables_relationships,fixed = TRUE))]," is nested within the combination of ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))], " and ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships) & !grepl(":",variables_relationships,fixed = TRUE))]))
      }
      else if (sum(!is.na(variables_relationships))==2 & sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==0){
        cat(paste0(mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1]," and ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2]," are crossed | ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][1]," and ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))][2]," are both nested within ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))]))
      }
      else if (sum(!is.na(variables_relationships))==1 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==1){
        cat(paste0(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1]," and ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2]," are crossed | ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))]," is nested within the combination of ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1]," and ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2]))
      }
      else if (sum(!is.na(variables_relationships))==1 &
               sum(grepl(":",variables_relationships[which(!is.na(variables_relationships))],fixed = TRUE))==0){
        cat(paste0(mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][1]," and ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships))][2]," are crossed | ",mlr3misc::extract_vars(form)$rhs[which(is.na(variables_relationships) & mlr3misc::extract_vars(form)$rhs!=variables_relationships[which(!is.na(variables_relationships))])]," and ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))]," are crossed | ",mlr3misc::extract_vars(form)$rhs[which(!is.na(variables_relationships))]," is nested within ",variables_relationships[which(!is.na(variables_relationships))]))
      }
    }
    cat("\n")
    }
  }
  # Just output the result table
  else {
    ifelse(print_ges==TRUE,ResultsToPrint<-(data.frame(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect!="(Intercept)",]) %>%
                                              dplyr::rename(Sig=.data$p..05,SS=.data$SSn,DF=.data$DFn,pval=.data$p) %>%
                                              dplyr::mutate(Sig=ifelse(.data$pval<sig_threshold,"*",""),
                                                            MSS=.data$SS/.data$DF) %>%
                                              base::rbind(data.frame(Effect="Residuals",SS=Results_of_ANOVA_test$SSd[1],SSd=NA,DF=Results_of_ANOVA_test$DFd[1],DFd=NA,Fval=NA,pval=NA,Sig="",ges=NA,MSS=unique(Results_of_ANOVA_test$SSd/Results_of_ANOVA_test$DFd))) %>%
                                              dplyr::mutate(SSd=NULL,DFd=NULL))[,c("Effect","SS","DF","MSS","Fval","pval","Sig","ges")],
           ResultsToPrint<-(data.frame(Results_of_ANOVA_test[Results_of_ANOVA_test$Effect!="(Intercept)",]) %>%
                              dplyr::rename(Sig=.data$p..05,SS=.data$SSn,DF=.data$DFn,pval=.data$p) %>%
                              dplyr::mutate(Sig=ifelse(.data$pval<sig_threshold,"*",""),
                                            MSS=.data$SS/.data$DF) %>%
                              base::rbind(data.frame(Effect="Residuals",SS=Results_of_ANOVA_test$SSd[1],SSd=NA,DF=Results_of_ANOVA_test$DFd[1],DFd=NA,Fval=NA,pval=NA,Sig="",ges=NA,MSS=unique(Results_of_ANOVA_test$SSd/Results_of_ANOVA_test$DFd))) %>%
                              dplyr::mutate(SSd=NULL,DFd=NULL))[,c("Effect","SS","DF","MSS","Fval","pval","Sig")])
    return(ResultsToPrint %>% `rownames<-`( NULL ))
  }
}
