################################################################################
#####                                Utils                                 #####
################################################################################

loss_mse <- function(h_hat, h)  mean((h_hat - h)^2)
loss_qlike <- function(h_hat, h) mean(log(h) + h_hat/h)
loss_qlike_ult <- function(h_hat, h) mean(h_hat/h - log(h_hat/h) - 1)
loss_mse_log <- function(h_hat, h) mean((log(h_hat) - log(h))^2)
loss_mse_sd <- function(h_hat, h) mean((sqrt(h_hat) - sqrt(h))^2)
loss_mse_prop <- function(h_hat, h) mean((h_hat / h - 1)^2)
loss_mae <- function(h_hat, h)  mean(abs(h_hat - h))
loss_mae_log <- function(h_hat, h)  mean(abs(log(h_hat) - log(h)))
loss_mae_sd <- function(h_hat, h) mean(abs(sqrt(h_hat) - sqrt(h)))
loss_mae_prop <- function(h_hat, h) mean(abs(h_hat / h - 1))


metrics <- function(h_hat, h) {
  c(loss_mse(h_hat, h),
    loss_qlike(h_hat, h),
    loss_qlike_ult(h_hat, h),
    loss_mse_log(h_hat, h),
    loss_mse_sd(h_hat, h),
    loss_mse_prop(h_hat, h),
    loss_mae(h_hat, h),
    loss_mae_log(h_hat, h),
    loss_mae_sd(h_hat, h),
    loss_mae_prop(h_hat, h))
}