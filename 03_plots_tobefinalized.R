plot_sensitivity(
    results,
    x = "Depth",
    title = NULL
)

plot_depth_test(test_depths)

plot_rarefaction_test(test_rare)

plot_absence_test(test_abs)



# matplot(
#     test_Nrare$Rarefactions,
#     test_Nrare[, c("Global","FHNY","LMDC")],
#     type = "b",
#     pch = 19,
#     lty = 1,
#     xlab = "Number of rarefactions",
#     ylab = "Number of core OTUs"
# )
# 
# legend(
#     "topright",
#     legend = c("Global","FHNY","LMDC"),
#     col = 1:3,
#     lty = 1,
#     pch = 19
# )
sample_richness <- colSums(stratified$otu > 0)
mean_richness <- tapply(
    sample_richness,
    stratified$meta$city,
    mean
)
sd_richness <- tapply(
    sample_richness,
    stratified$meta$city,
    sd
)
data.frame(
    City = names(mean_richness),
    Mean = round(mean_richness, 1),
    SD = round(sd_richness, 1),
    Total_OTUs = total_detected[names(mean_richness)]
)

# test_abs$FHNY_per_sample <- test_abs$FHNY / mean_richness["FHNY"]
# test_abs$LMDC_per_sample <- test_abs$LMDC / mean_richness["LMDC"]
# 
# matplot(
#     test_abs$MaxAbs,
#     test_abs[,c("Global","FHNY","LMDC")],
#     type="b",
#     pch=19,
#     lty=1,
#     xlab="Allowed absences",
#     ylab="Number of core OTUs"
# )
# 
# legend(
#     "topleft",
#     legend=c("Global","FHNY","LMDC"),
#     col=1:3,
#     pch=19,
#     lty=1
# )
# text(
#     x = 5,
#     y = test_abs$FHNY[nrow(test_abs)],
#     labels = sprintf("%.2f%%", test_abs$FHNY_per_sample[nrow(test_abs)]),
#     col = "red",
#     pos = 4, xpd = NA
# )
# 
# text(
#     x = 5,
#     y = test_abs$LMDC[nrow(test_abs)],
#     labels = sprintf("%.2f%%", test_abs$LMDC_per_sample[nrow(test_abs)]),
#     col = "green4",
#     pos = 4, xpd = NA
# )
