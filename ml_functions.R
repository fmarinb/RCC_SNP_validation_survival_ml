library(caret)
library(MLeval)
library(tidyverse)
library(pROC)


ml_inference <- function(
  datalist,
  y,
  ctrl_method,
  number = NULL,
  metric = "ROC",
  pos_class = NULL,
  tuneLength = 10
) {
  # helpers -------------------------
  get_train_cm <- function(fit, positive) {
    p <- fit$pred

    # si hay hiperparámetros, filtrar al mejor
    if (!is.null(fit$bestTune)) {
      for (nm in names(fit$bestTune)) {
        p <- p[p[[nm]] == fit$bestTune[[nm]], ]
      }
    }

    caret::confusionMatrix(
      data = factor(p$pred, levels = levels(p$obs)),
      reference = p$obs,
      positive = positive
    )
  }

  f1_from_cm <- function(cm) {
    ppv <- cm$byClass["Pos Pred Value"]
    rec <- cm$byClass["Sensitivity"]
    ifelse((ppv + rec) == 0, 0, 2 * ppv * rec / (ppv + rec))
  }
  # ---------------------------------

  results_list <- lapply(datalist, function(data) {
    set.seed(123)

    data <- dplyr::select(data, -dplyr::any_of(c("ID", "class")))

    trainIndex <- caret::createDataPartition(data[[y]], p = 0.70, list = FALSE)
    train_data <- data[trainIndex, ]
    test_data <- data[-trainIndex, ]

    fml_y <- stats::as.formula(sprintf("%s ~ .", y))

    ctrl <- caret::trainControl(
      method = ctrl_method,
      number = number,
      summaryFunction = caret::twoClassSummary,
      classProbs = TRUE,
      savePredictions = "final"
    )

    message("start Random Forests")
    rf <- caret::train(
      fml_y,
      data = train_data,
      method = "rf",
      trControl = ctrl,
      metric = metric,
      tuneLength = tuneLength
    )

    message("Random Forests finished: start C5.0 tree")
    c5.0 <- caret::train(
      fml_y,
      data = train_data,
      method = "C5.0",
      trControl = ctrl,
      metric = metric,
      tuneLength = tuneLength
    )

    message("C5.0 tree finished: start glm Elastic Net")
    glmnet <- caret::train(
      fml_y,
      data = train_data,
      method = "glmnet",
      trControl = ctrl,
      metric = metric,
      tuneLength = tuneLength
    )

    message("glm Elastic Net finished: start radial kernel SVM")
    svm <- caret::train(
      fml_y,
      data = train_data,
      method = "svmRadial",
      trControl = ctrl,
      metric = metric,
      tuneLength = tuneLength
    )

    message("SVM finished: start logistic regression")
    logitglm <- caret::train(
      fml_y,
      data = train_data,
      method = "glm",
      family = "binomial",
      trControl = ctrl,
      metric = metric
    )

    # ===== Evaluación en TEST (AUC/ROC) =====
    yy <- test_data[[y]]
    if (!is.factor(yy)) {
      yy <- factor(yy)
    }
    if (is.null(pos_class)) {
      pos_class <- levels(yy)[2]
    }

    models <- list(
      "Logistic regression" = logitglm,
      "Elastic net" = glmnet,
      "C5.0 Decision tree" = c5.0,
      "Random Forests" = rf,
      "Radial kernel SVM" = svm
    )

    pred_prob <- function(fit) {
      predict(fit, newdata = test_data, type = "prob")[, pos_class]
    }
    probs_test <- lapply(models, pred_prob)

    roc_list <- lapply(probs_test, function(p) {
      pROC::roc(response = yy, predictor = p, quiet = TRUE)
    })
    auc_vec <- vapply(roc_list, pROC::auc, numeric(1))

    # ===== Matrices de confusión y F1 =====
    # TRAIN (OOF del resampling) – no altera nada entrenado
    cm_train <- lapply(models, get_train_cm, positive = pos_class)
    f1_train <- vapply(cm_train, f1_from_cm, numeric(1))

    # TEST
    pred_class_test <- lapply(
      models,
      function(fit) {
        predict(fit, newdata = test_data)
      } # clases
    )
    cm_test <- mapply(
      function(pred) {
        caret::confusionMatrix(
          data = pred,
          reference = yy,
          positive = pos_class
        )
      },
      pred_class_test,
      SIMPLIFY = FALSE
    )
    f1_test <- vapply(cm_test, f1_from_cm, numeric(1))

    # ===== OUTPUT =====
    list(
      models = models,
      auc_test = auc_vec,
      probs_test = probs_test,
      roc_list = roc_list,
      cm_train = cm_train,
      f1_train = f1_train,
      cm_test = cm_test,
      f1_test = f1_test,
      train_data = train_data,
      test_data = test_data
    )
  })

  return(results_list)
}


test_results <- function(ml_model_list, ml_test_df, y) {
  test_pred <- function(model) {
    p_prob <- predict(model, newdata = ml_test_df, type = 'prob')
    p_class <- predict(model, newdata = ml_test_df, type = 'raw')

    cm <- confusionMatrix(p_class, ml_test_df[[y]])

    preds <- data.frame(
      Healthy = p_prob[[1]],
      Tumoral = p_prob[[2]],
      obs = ml_test_df[[y]]
    )

    roc_obj <- roc(ml_test_df[[y]], preds$Healthy)
    auc_value <- auc(roc_obj)

    return(list(
      confusion_matrix = cm,
      predictions = preds,
      auc = auc_value
    ))
  }

  out <- lapply(ml_model_list, test_pred)
  names(out) <- names(ml_model_list)

  return(out)
}
