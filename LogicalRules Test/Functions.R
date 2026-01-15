library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

prepare_density_data <- function(data_list,
                                 factors = list(),
                                 subject_col = "subjects",
                                 correct_col = "corr",
                                 rt_col = "rt",
                                 group_level=FALSE,
                                 split_correct = TRUE,
                                 design_label = NULL) {
  # Apply factor filters to each data frame in the list
  data_list <- map(data_list, function(df) {
    for (factor_name in names(factors)) {
      factor_val <- factors[[factor_name]]
      df <- df[df[[factor_name]] == factor_val, , drop = FALSE]
    }
    df
  })
  # All subjects present in all sources
  subj_list <- map(data_list, ~unique(.x[[subject_col]]))
  subjects <- factor(Reduce(intersect, subj_list))
  if(!group_level) {
    results <- map_dfr(subjects, function(subj) {
      # For each source/model
      map_dfr(names(data_list), function(source_name) {
        data_subj <- data_list[[source_name]][data_list[[source_name]][[subject_col]] == subj, ]
        if (nrow(data_subj) == 0) return(NULL)
        
        if (split_correct) {
          dq_c <- data_subj[[rt_col]][data_subj[[correct_col]] == 1]
          dq_e <- data_subj[[rt_col]][data_subj[[correct_col]] == 0]
          
          bind_rows(
            tibble(
              Subject = subj,
              RT = as.numeric(dq_c),
              Source = source_name,
              Subset = "correct"
            ),
            tibble(
              Subject = subj,
              RT = as.numeric(dq_e),
              Source = source_name,
              Subset = "error"
            )
          )
        } else {
          dq <- data_subj[[rt_col]]
          tibble(
            Subject = subj,
            RT = as.numeric(dq),
            Source = source_name,
            Subset = "all"
          )
        }
      })
    }) 
  } else {
    results = map_dfr(names(data_list), function(source_name) {
      data <- data_list[[source_name]]
      if (nrow(data) == 0) return(NULL)
      if (nrow(data) == 0) return(NULL)
      
      if (split_correct) {
        dq_c <- data[[rt_col]][data[[correct_col]] == 1]
        dq_e <- data[[rt_col]][data[[correct_col]] == 0]
        
        bind_rows(
          tibble(
            RT = as.numeric(dq_c),
            Source = source_name,
            Subset = "correct"
          ),
          tibble(
            RT = as.numeric(dq_e),
            Source = source_name,
            Subset = "error"
          )
        )
      } else {
        dq <- data[[rt_col]]
        tibble(
          RT = as.numeric(dq),
          Source = source_name,
          Subset = "all"
        )
      }
    })
  }
  
  if (!is.null(design_label)) {
    results <- dplyr::mutate(results, DesignCell = design_label)
  }
  
  return(results)
}

compute_quantile_points_by_subject <- function(data_list,
                                               subject_col = "subjects",
                                               group_level=FALSE,
                                               factors = list(),
                                               correct_col = "corr",
                                               rt_col = "rt",
                                               quantiles = seq(0.05, 0.95, 0.05),
                                               split_correct = TRUE,
                                               postn=NULL,
                                               design_label = NULL) {
  # Apply factor filters to each data frame in the list
  data_list <- map(data_list, function(df) {
    for (factor_name in names(factors)) {
      factor_val <- factors[[factor_name]]
      df <- df[df[[factor_name]] == factor_val, , drop = FALSE]
    }
    df
  })
  # All subjects present in all sources
  subj_list <- map(data_list, ~unique(.x[[subject_col]]))
  subjects <- Reduce(intersect, subj_list)
  if(!group_level) {
    results <- map_dfr(subjects, function(subj) {
      # For each source/model
      map_dfr(names(data_list), function(source_name) {
        data_subj <- data_list[[source_name]][data_list[[source_name]][[subject_col]] == subj, ]
        if (nrow(data_subj) == 0) return(NULL)
        
        # Proportion correct (can override with p_correct_list)
        pc =mean(data_subj[[correct_col]], na.rm = TRUE)
        
        if (split_correct) {
          dq_c <- quantile(data_subj[[rt_col]][data_subj[[correct_col]] == 1], quantiles, na.rm = TRUE)
          dq_e <- quantile(data_subj[[rt_col]][data_subj[[correct_col]] == 0], quantiles, na.rm = TRUE)
          
          bind_rows(
            tibble(
              Subject = subj,
              Quantile = quantiles,
              RT = as.numeric(dq_c),
              Source = source_name,
              Subset = "correct",
              RT_weighted = as.numeric(dq_c) * pc
            ),
            tibble(
              Subject = subj,
              Quantile = quantiles,
              RT = as.numeric(dq_e),
              Source = source_name,
              Subset = "error",
              RT_weighted = as.numeric(dq_e) * (1 - pc)
            )
          )
        } else {
          dq <- quantile(data_subj[[rt_col]], quantiles, na.rm = TRUE)
          tibble(
            Subject = subj,
            Quantile = quantiles,
            RT = as.numeric(dq),
            Source = source_name,
            Subset = "all",
            RT_weighted = as.numeric(dq) # weighting logic can be updated
          )
        }
      })
    }) 
  } else {
    results = map_dfr(names(data_list), function(source_name) {
      data <- data_list[[source_name]]
      if (nrow(data) == 0) return(NULL)
      if (nrow(data) == 0) return(NULL)
      
      # Proportion correct (can override with p_correct_list)
      pc = mean(data[[correct_col]], na.rm = TRUE)
      
      if (split_correct) {
        dq_c <- quantile(data[[rt_col]][data[[correct_col]] == 1], quantiles, na.rm = TRUE)
        dq_e <- quantile(data[[rt_col]][data[[correct_col]] == 0], quantiles, na.rm = TRUE)
        
        bind_rows(
          tibble(
            Quantile = quantiles,
            RT = as.numeric(dq_c),
            Source = source_name,
            Subset = "correct",
            RT_weighted = as.numeric(dq_c) * pc
          ),
          tibble(
            Quantile = quantiles,
            RT = as.numeric(dq_e),
            Source = source_name,
            Subset = "error",
            RT_weighted = as.numeric(dq_e) * (1 - pc)
          )
        )
      } else {
        dq <- quantile(data[[rt_col]], quantiles, na.rm = TRUE)
        tibble(
          Quantile = quantiles,
          RT = as.numeric(dq),
          Source = source_name,
          Subset = "all",
          RT_weighted = as.numeric(dq) # weighting logic can be updated
        )
      }
    })
  }
  
  if (!is.null(design_label)) {
    results <- dplyr::mutate(results, DesignCell = design_label)
  }
  
  return(results)
}

compute_accuracy_by_subject <- function(
    data_list,
    subject_col = "subjects",
    group_level = FALSE,
    factors = list(),
    correct_col = "corr",
    design_label = NULL,
    p_correct_list = NULL  # Optionally, a named vector/list to override computed means
) {
  # Apply factor filters to each data frame in the list
  data_list <- purrr::map(data_list, function(df) {
    for (factor_name in names(factors)) {
      factor_val <- factors[[factor_name]]
      df <- df[df[[factor_name]] == factor_val, , drop = FALSE]
    }
    df
  })
  
  # Get subject intersection, if subject-level
  if (!group_level) {
    subj_list <- map(data_list, ~unique(.x[[subject_col]]))
    subjects <- Reduce(intersect, subj_list)
    results <- map_dfr(subjects, function(subj) {
      map_dfr(names(data_list), function(source_name) {
        df_subj <- data_list[[source_name]][data_list[[source_name]][[subject_col]] == subj, , drop = FALSE]
        if (nrow(df_subj) == 0) return(NULL)
        pc <- if (!is.null(p_correct_list) && !is.null(p_correct_list[[source_name]])) {
          p_correct_list[[source_name]]
        } else {
          mean(df_subj[[correct_col]], na.rm = TRUE)
        }
        tibble(
          Subject = subj,
          Accuracy = as.numeric(pc),
          Source = source_name,
          Subset = "all"
        )
      })
    })
  } else {
    results <- map_dfr(names(data_list), function(source_name) {
      df <- data_list[[source_name]]
      pc <- if (!is.null(p_correct_list) && !is.null(p_correct_list[[source_name]])) {
        p_correct_list[[source_name]]
      } else {
        mean(df[[correct_col]], na.rm = TRUE)
      }
      tibble(
        Accuracy = as.numeric(pc),
        Source = source_name,
        Subset = "all"
      )
    })
  }
  
  if (!is.null(design_label)) {
    results <- dplyr::mutate(results, DesignCell = design_label)
  }
  
  return(results)
}
