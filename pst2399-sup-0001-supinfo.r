###################
## TGI tutorial  ##
###################

#########################################################################################################################################################################
## Box 1:
## A simulated dataset was created with the following characteristics:
## One vehicle group, 2 single agent therapies (A and B) group, and one combination therapy A+B group ('group.name'). Each group starts with 15 animals ('n.animals').
## Day 0 is when animals are assigned to groups and treatment starts, if any. The study lasts for 18 days ('duration').
## All animals start with a tumor volume of 100 mm^3 ('initial.tumor.volume') and have their tumor volumes measured every 3 days. The measuring error amounts to CV = 25% ('residual.cv').
## The average tumor growth rates (in %) for the groups are: Vehicle = 5; Therapy_A = 4.5; Therapy_B = 3.5; Therapy_AB = 1 ('group.growth.rate.day.perc') with intra-group CV=10% ('growth.rate.cv')
## Animals whose tumor volume reach 4000 mm^3 are removed from subsequent days of the study, if any ('max.tumor.volume').

simulate.exp.growth <- function(group.name, duration, initial.tumor.volume, group.growth.rate.day.perc, growth.rate.cv, residual.cv, n.animals, max.tumor.volume) {
    group <- list()
    for (i in 1:n.animals) {
        ln.tumor.volume <- numeric()
        ln.tumor.volume[1] <- rnorm(n = 1, mean = log(initial.tumor.volume), sd = sqrt(log((residual.cv^2) + 1)))
        days <- seq(0, duration, 3)
        growth.rate.individual <- rnorm(n = (length(days) - 1), mean = group.growth.rate.day.perc / 100, sd = group.growth.rate.day.perc * growth.rate.cv / 100)

        for (t in 2:length(days)) {
            ln.tumor.volume[t] <- rnorm(n = 1, mean = ln.tumor.volume[t - 1] + (growth.rate.individual[t - 1] * days[t]), sd = sqrt(log((residual.cv^2) + 1)))
        }
        welfare.missing <- which(exp(ln.tumor.volume) > max.tumor.volume)
        if (length(welfare.missing) > 1) ln.tumor.volume[min(welfare.missing[-1]):length(ln.tumor.volume)] <- NA
        group[[i]] <- data.frame(Day = days, Animal = paste(group.name, i, sep = "_"), Group = group.name, TV = exp(ln.tumor.volume))
    }
    group <- plyr::rbind.fill(group)
    group
}

set.seed(92121)
vehicle.group <- simulate.exp.growth(
    group.name = "Control", duration = 18, initial.tumor.volume = 100, group.growth.rate.day.perc = 5,
    growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 15, max.tumor.volume = 4000
)
therapyA.group <- simulate.exp.growth(
    group.name = "Therapy_A", duration = 18, initial.tumor.volume = 100, group.growth.rate.day.perc = 4,
    growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 15, max.tumor.volume = 4000
)
therapyB.group <- simulate.exp.growth(
    group.name = "Therapy_B", duration = 18, initial.tumor.volume = 100, group.growth.rate.day.perc = 3,
    growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 15, max.tumor.volume = 4000
)
therapyAB.group <- simulate.exp.growth(
    group.name = "Therapy_Combo_AB", duration = 18, initial.tumor.volume = 100, group.growth.rate.day.perc = 1.3,
    growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 15, max.tumor.volume = 4000
)
TGI.dataset <- rbind(vehicle.group, therapyA.group, therapyB.group, therapyAB.group)

head(TGI.dataset)

##   Day    Animal   Group       TV
## 1   0 Control_1 Control 104.8260
## 2   3 Control_1 Control 125.1373
## 3   6 Control_1 Control 155.8472
## 4   9 Control_1 Control 238.6878
## 5  12 Control_1 Control 335.1191
## 6  15 Control_1 Control 729.3432



#########################################################################################################################################################################
## Box 1:
## Function to obtain tumor growth summary statistics:
##
## summary.stats.TGI()
## Required argument is a TGI dataset formatted as 'TGI.dataset'
## Outputs geometric means and respective standard errors and 95% confidence intervals per animal group per day.

summary.stats.TGI <- function(dataset) {
    dataset.list <- split(dataset, list(dataset$Group, dataset$Day), drop = TRUE)
    dataset.list <- lapply(dataset.list, function(x) {
        new.x <- as.data.frame(unique(x[, c("Day", "Group")]))
        new.x$est <- exp(mean(log(x$TV), na.rm = TRUE))
        new.x$se <- sqrt(var(log(x$TV), na.rm = TRUE) / nrow(x[!is.na(x$TV), ]))
        new.x$lb95.est <- exp(log(new.x$est) + qnorm(0.025) * new.x$se)
        new.x$ub95.est <- exp(log(new.x$est) + qnorm(0.975) * new.x$se)
        new.x
    })
    dataset.list <- plyr::rbind.fill(dataset.list)
    dataset.list
}

#########################################################################################################################################################################
## Box 1:
## Function to plot tumor growth curves:
##
## plot.TGI.dataset()
## Required argumes is a data-frame containing averages and confidence intervals in the same format as output by 'summary.stats.TGI' function
## Depicts a plot with tumor growth curves in the original scale (default) or in natural ln-scale (if ln.scale = "yes")
## Tumor growth curves are plotted in colors (default) or in gray-shades (if grey.shades = "yes")

plot.TGI.dataset <- function(summ.stats, ln.scale = "no", grey.shades = "no") {
    summ.stats$Group <- as.factor(summ.stats$Group)
    summ.stats$group.color <- viridis::viridis(nlevels(summ.stats$Group) + 1)[as.numeric(summ.stats$Group)] ## requires package 'viridis' to be installed
    if (grey.shades == "yes") summ.stats$group.color <- grDevices::gray.colors(n = nlevels(summ.stats$Group) + 1)[as.numeric(summ.stats$Group)] ## requires package 'grDevices' to be installed
    summ.stats$group.pch <- c(15, 16, 17, 18)[as.numeric(summ.stats$Group)]
    y.ticks <- seq(0, max(summ.stats$ub95.est, na.rm = TRUE), 500)

    if (ln.scale == "yes") {
        summ.stats$est <- log(summ.stats$est)
        summ.stats$lb95.est <- log(summ.stats$lb95.est)
        summ.stats$ub95.est <- log(summ.stats$ub95.est)
        y.ticks <- log(y.ticks)
    }
    par(mar = c(5, 6, 3, 1))
    plotrix::plotCI(x = summ.stats$Day, y = summ.stats$est, li = summ.stats$lb95.est, ui = summ.stats$ub95.est, xaxt = "n", xlab = "Days", yaxt = "n", ylab = "", cex.lab = 1.3, font.lab = 2, main = "Tumor Growth Curves", cex.main = 2.5, font.main = 2, pch = NA, slty = 0, lwd = 1.7)
    axis(1, at = unique(summ.stats$Day), label = unique(summ.stats$Day))
    abline(v = unique(summ.stats$Day), lty = 3, col = "gray80")
    if (ln.scale == "yes") {
        axis(2, at = y.ticks, label = format(exp(y.ticks), sci = T), las = 2)
    } else {
        axis(2, at = y.ticks, label = format(y.ticks, sci = T), las = 2)
    }
    abline(h = y.ticks, lty = 3, col = "gray80")
    mtext(side = 2, text = "Tumor Volume (mm^3)", line = 5, font = 2, cex = 1.3)
    group.data <- split(summ.stats, summ.stats$Group)
    invisible(lapply(group.data, function(x) {
        plotrix::plotCI(x = x$Day, y = x$est, li = x$lb95.est, ui = x$ub95.est, pch = unique(x$group.pch), col = unique(x$group.color), add = TRUE, cex = 1.5)
        lines(x = x$Day, y = x$est, col = unique(x$group.color), lwd = 2)
    }))
    legend.attributes <- unique(summ.stats[, c("Group", "group.color", "group.pch")])
    legend("topleft", legend = legend.attributes$Group, col = legend.attributes$group.color, pch = legend.attributes$group.pch, lwd = 1.4, cex = 1.3, bty = "n")
}

TGI.stats <- summary.stats.TGI(dataset = TGI.dataset)
plot.TGI.dataset(TGI.stats)
plot.TGI.dataset(TGI.stats, ln.scale = "yes")

#########################################################################################################################################################################
## Box 2:
## Calculate TGI from summary statistics:
##
## calculate.TGI.from.data()
## Requires a dataset in the same format as TGI.dataset ('dataset'). TGI or TGI-delta can be calculate ('type'). Argument accepts one or both.
## Time-points for which the TGI will be calculated should be passed to 'time.points'. If TGI-delta is calculated, 'baseline' needs a valid baseline day
## The names of the treatment ('treatment.name') and control ('vehicle.name') groups are required.

calculate.TGI.from.data <- function(dataset, type = c("TGI", "TGI.delta"), time.points = NULL, baseline = NULL, treatment.name, vehicle.name) {
    if (any(type %in% "TGI.delta") & is.null(baseline)) stop("Please assign a baseline time point for TGI delta calculations")
    vehicle.data <- subset(dataset, Group %in% vehicle.name)
    if (!is.null(time.points)) vehicle.data <- subset(vehicle.data, Day %in% c(baseline, time.points))
    if (!is.null(baseline)) {
        vehicle.baseline <- subset(vehicle.data, Day %in% baseline)
        vehicle.data <- subset(vehicle.data, !(Day %in% baseline))
    }
    treatment.data <- subset(dataset, Group %in% treatment.name)
    if (!is.null(time.points)) treatment.data <- subset(treatment.data, Day %in% c(baseline, time.points))
    if (!is.null(baseline)) {
        treatment.baseline <- subset(treatment.data, Day %in% baseline)
        treatment.data <- subset(treatment.data, !(Day %in% baseline))
    }

    vehicle.means <- as.data.frame(aggregate(x = vehicle.data[, c("Day", "TV")], by = list(vehicle.data$Day), function(x) mean(x, na.rm = TRUE))[, -1])
    names(vehicle.means)[names(vehicle.means) == "TV"] <- "TV.vehicle"
    if (!is.null(baseline)) vehicle.baseline.means <- as.data.frame(aggregate(x = vehicle.baseline[, c("Day", "TV")], by = list(vehicle.baseline$Day), function(x) mean(x, na.rm = TRUE))[, -1])

    treatment.means <- as.data.frame(aggregate(x = treatment.data[, c("Day", "TV")], by = list(treatment.data$Day), function(x) mean(x, na.rm = TRUE))[, -1])
    names(treatment.means)[names(treatment.means) == "TV"] <- "TV.treatment"
    if (!is.null(baseline)) treatment.baseline.means <- as.data.frame(aggregate(x = treatment.baseline[, c("Day", "TV")], by = list(treatment.baseline$Day), function(x) mean(x, na.rm = TRUE))[, -1])

    TGI <- NULL
    TGI.delta <- NULL

    if (any(type %in% "TGI")) {
        TGI.data <- merge(vehicle.means, treatment.means, by = "Day", sort = FALSE)
        TGI.data$TGI <- (1 - (TGI.data$TV.treatment / TGI.data$TV.vehicle)) * 100
        TGI.data$Contrast <- paste(treatment.name, "vs.", vehicle.name)
        TGI <- TGI.data[, c("Contrast", "Day", "TGI")]
    }

    if (any(type %in% "TGI.delta")) {
        vehicle.means$TV.vehicle.adjusted <- vehicle.means$TV.vehicle - vehicle.baseline.means$TV
        treatment.means$TV.treatment.adjusted <- treatment.means$TV.treatment - treatment.baseline.means$TV
        TGI.delta.data <- merge(vehicle.means, treatment.means, by = "Day", sort = FALSE)
        TGI.delta.data$TGI.delta <- (1 - (TGI.delta.data$TV.treatment.adjusted / TGI.delta.data$TV.vehicle.adjusted)) * 100
        TGI.delta.data$Contrast <- paste(treatment.name, "vs.", vehicle.name)
        TGI.delta <- TGI.delta.data[, c("Contrast", "Day", "TGI.delta")]
    }

    return(list(TGI = TGI, TGI.delta = TGI.delta))
}

TGI.therapyB.vs.vehicle <- calculate.TGI.from.data(dataset = TGI.dataset, type = "TGI", time.points = c(18, 21, 24, 27), treatment.name = "Therapy_B", vehicle.name = "Control")
TGI.delta.therapyB.vs.vehicle <- calculate.TGI.from.data(dataset = TGI.dataset, type = "TGI.delta", time.points = c(18, 21, 24, 27), baseline = 0, treatment.name = "Therapy_B", vehicle.name = "Control")

#########################################################################################################################################################################
## Box 3:
## Analysis of continuous data at single time-points versus vehicle (one-sided) or versus another treatment (two-sided):
##
## TGI.single.time()
## Requires a dataset in the same format as TGI.dataset ('dataset'). Estimates TGI for a single time-point ('time.point').
## TGI estimates can be either adjusted or unadjusted ('baseline.adjusted') for baseline day ('baseline').
## The names of the treatments ('treatments.numerator') being compared to the ('treatment.denominator') groups are required.
## One-sided or two-sided tests ('p.tails')

TGI.single.time <- function(dataset, time.point, baseline.adjusted = "no", baseline = NULL, treatments.numerator, treatment.denominator, p.tails = 1) {
    if (baseline.adjusted == "yes" & is.null(baseline)) stop("Please either assign a baseline time point for TGI calculations or set 'baseline.adjusted' to 'no'.")
    TV.data <- subset(dataset, Day %in% c(baseline, time.point) & Group %in% c(treatments.numerator, treatment.denominator))
    TV.data$lnTV <- log(TV.data$TV)
    if (baseline.adjusted == "yes") {
        TV.baseline <- subset(TV.data, Day %in% baseline)[, c("Animal", "lnTV")]
        names(TV.baseline)[names(TV.baseline) == "lnTV"] <- "lnTV.baseline"
        TV.data <- subset(TV.data, !(Day %in% baseline))
        TV.data <- merge(TV.data, TV.baseline, by = "Animal", sort = FALSE)
    }
    TV.data$Group <- factor(TV.data$Group, levels = c(treatment.denominator, treatments.numerator))

    if (baseline.adjusted == "no") model.fit <- lm(lnTV ~ Group, data = TV.data)
    if (baseline.adjusted == "yes") model.fit <- lm(lnTV ~ lnTV.baseline + Group, data = TV.data)

    model.estimates <- as.data.frame(coef(summary(model.fit)))
    model.estimates <- model.estimates[paste("Group", treatments.numerator, sep = ""), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
    model.estimates$TGI <- (1 - exp(model.estimates$Estimate)) * 100
    if (p.tails == 1) {
        model.estimates$lb95.TGI <- NA
        model.estimates$ub95.TGI <- (1 - exp(model.estimates$Estimate + qt(0.95, model.fit$df.residual) * model.estimates$"Std. Error")) * 100
        model.estimates$pvalue <- pt(model.estimates$"t value", model.fit$df.residual)
    }
    if (p.tails == 2) {
        model.estimates$lb95.TGI <- (1 - exp(model.estimates$Estimate + qt(0.025, model.fit$df.residual) * model.estimates$"Std. Error")) * 100
        model.estimates$ub95.TGI <- (1 - exp(model.estimates$Estimate + qt(0.975, model.fit$df.residual) * model.estimates$"Std. Error")) * 100
        model.estimates$pvalue <- 2 * (1 - pt(abs(model.estimates$"t value"), model.fit$df.residual))
    }
    model.estimates$Contrast <- paste(treatments.numerator, "vs.", treatment.denominator)
    model.estimates$Day <- time.point
    model.estimates$Df <- model.fit$df.residual
    row.names(model.estimates) <- NULL
    model.estimates <- model.estimates[, c("Contrast", "Day", "TGI", "lb95.TGI", "ub95.TGI", "pvalue", "Df")]

    return(model.estimates)
}

TGI.model.therapyA.or.B.vs.vehicle.day18 <- TGI.single.time(dataset = TGI.dataset, time.point = 18, baseline.adjusted = "no", treatments.numerator = c("Therapy_A", "Therapy_B"), treatment.denominator = "Control")
TGI.baseline.adjusted.model.therapyA.or.B.vs.vehicle.day18 <- TGI.single.time(dataset = TGI.dataset, time.point = 18, baseline.adjusted = "yes", baseline = 0, treatments.numerator = c("Therapy_A", "Therapy_B"), treatment.denominator = "Control")

TGI.model.therapyA.vs.therapyB.day18 <- TGI.single.time(dataset = TGI.dataset, time.point = 18, baseline.adjusted = "no", treatments.numerator = "Therapy_A", treatment.denominator = "Therapy_B", p.tails = 2)
TGI.baseline.adjusted.model.therapyA.vs.therapyB.day18 <- TGI.single.time(dataset = TGI.dataset, time.point = 18, baseline.adjusted = "yes", baseline = 0, treatments.numerator = "Therapy_A", treatment.denominator = "Therapy_B", p.tails = 2)

#########################################################################################################################################################################
## Box 4:
## Analysis of continuous data at multiple time-points versus vehicle (one-sided) or versus another treatment (two-sided):
##
## TGI.multiple.time()
## Requires a dataset in the same format as TGI.dataset ('dataset'). Estimates TGI for multiple time-points ('time.points').
## TGI estimates can be either adjusted or unadjusted ('baseline.adjusted') for baseline day ('baseline').
## The names of the treatments ('treatments.numerator') being compared to the ('treatment.denominator') groups are required.
## One-sided or two-sided tests ('p.tails')

TGI.multiple.time <- function(dataset, time.points, baseline.adjusted = "no", baseline = NULL, treatments.numerator, treatment.denominator, p.tails = 1) {
    if (length(time.points) < 2) stop("Please use 'TGI.single.time' function.")
    if (baseline.adjusted == "yes" & is.null(baseline)) stop("Please either assign a baseline time point for TGI calculations or set 'baseline.adjusted' to 'no'.")
    TV.data <- subset(dataset, Day %in% c(baseline, time.points) & Group %in% c(treatments.numerator, treatment.denominator))
    TV.data$lnTV <- log(TV.data$TV)
    if (baseline.adjusted == "yes") {
        TV.baseline <- subset(TV.data, Day %in% baseline)[, c("Animal", "lnTV")]
        names(TV.baseline)[names(TV.baseline) == "lnTV"] <- "lnTV.baseline"
        TV.data <- subset(TV.data, !(Day %in% baseline))
        TV.data <- merge(TV.data, TV.baseline, by = "Animal", sort = FALSE)
    }
    TV.data$Group <- factor(TV.data$Group, levels = c(treatment.denominator, treatments.numerator))
    TV.data$Day <- factor(TV.data$Day, levels = sort(unique(TV.data$Day)))
    TV.data.list <- split(TV.data, TV.data$Day, drop = TRUE)

    model.estimates <- lapply(TV.data.list, function(x) {
        if (baseline.adjusted == "no") model.fit <- lm(lnTV ~ Group, data = x)
        if (baseline.adjusted == "yes") model.fit <- lm(lnTV ~ lnTV.baseline + Group, data = x)

        model.output <- as.data.frame(coef(summary(model.fit)))
        model.output <- model.output[paste("Group", treatments.numerator, sep = ""), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
        model.output$TGI <- (1 - exp(model.output$Estimate)) * 100
        if (p.tails == 1) {
            model.output$lb95.TGI <- NA
            model.output$ub95.TGI <- (1 - exp(model.output$Estimate + qt(0.95, model.fit$df.residual) * model.output$"Std. Error")) * 100
            model.output$pvalue <- pt(model.output$"t value", model.fit$df.residual)
        }
        if (p.tails == 2) {
            model.output$lb95.TGI <- (1 - exp(model.output$Estimate + qt(0.025, model.fit$df.residual) * model.output$"Std. Error")) * 100
            model.output$ub95.TGI <- (1 - exp(model.output$Estimate + qt(0.975, model.fit$df.residual) * model.output$"Std. Error")) * 100
            model.output$pvalue <- 2 * (1 - pt(abs(model.output$"t value"), model.fit$df.residual))
        }
        model.output$Day <- as.character(as.vector(unique(x$Day)))
        model.output$Df <- model.fit$df.residual
        model.output
    })

    model.estimates <- plyr::rbind.fill(model.estimates)
    model.estimates$Contrast <- paste(treatments.numerator, "vs.", treatment.denominator)
    model.estimates <- model.estimates[, c("Contrast", "Day", "TGI", "lb95.TGI", "ub95.TGI", "pvalue", "Df")]

    return(model.estimates)
}

TGI.model.all.therapies.vs.vehicle.all.days <- TGI.multiple.time(dataset = TGI.dataset, time.points = c(0, 3, 6, 9, 12, 15, 18), baseline.adjusted = "no", treatments.numerator = c("Therapy_A", "Therapy_B", "Therapy_Combo_AB"), treatment.denominator = "Control")
TGI.baseline.adjusted.model.all.therapies.vs.vehicle.all.days <- TGI.multiple.time(dataset = TGI.dataset, time.points = c(3, 6, 9, 12, 15, 18), baseline.adjusted = "yes", baseline = 0, treatments.numerator = c("Therapy_A", "Therapy_B", "Therapy_Combo_AB"), treatment.denominator = "Control")
TGI.baseline.adjusted.model.therapyA.vs.therapyB.all.days <- TGI.multiple.time(dataset = TGI.dataset, time.points = c(3, 6, 9, 12, 15, 18), baseline.adjusted = "yes", baseline = 0, treatments.numerator = "Therapy_A", treatment.denominator = "Therapy_B", p.tails = 2)


#########################################################################################################################################################################
## Box 5:
## Analysis of growth curves exponential:
##
## Tumor.growth.curves.exponential()
## Requires a dataset in the same format as TGI.dataset ('dataset').
## Compares tumor growth rates "TGR" of treatments ('treatments.numerator') versus control or other treatments ('treatment.denominator').
## One-sided or two-sided tests ('p.tails')

Tumor.growth.curves.exponential <- function(dataset, treatments.numerator, treatment.denominator, p.tails = 1) {
    TV.data <- subset(dataset, Group %in% c(treatments.numerator, treatment.denominator))
    TV.data$lnTV <- log(TV.data$TV)
    TV.data$Group <- factor(TV.data$Group, levels = c(treatment.denominator, treatments.numerator))
    TV.data$Day <- as.numeric(TV.data$Day)
    model.fit <- nlme::lme(lnTV ~ Group * Day, random = ~ 1 + Day | Animal, data = TV.data, na.action = na.omit, method = "REML")
    model.estimates <- as.data.frame(summary(model.fit)$tTable)
    model.estimates <- model.estimates[grep(":Day", rownames(model.estimates)), ]
    names(model.estimates)[names(model.estimates) == "Value"] <- "TGR"
    model.estimates$lb95.TGR <- NA
    model.estimates$ub95.TGR <- NA
    model.estimates$pvalue <- NA
    for (i in 1:nrow(model.estimates)) {
        if (p.tails == 1) {
            model.estimates$ub95.TGR[i] <- model.estimates$TGR[i] + qt(0.95, model.estimates$DF[i]) * model.estimates$"Std.Error"[i]
            model.estimates$pvalue[i] <- pt(model.estimates$"t-value"[i], model.estimates$DF[i])
        }
        if (p.tails == 2) {
            model.estimates$lb95.TGR[i] <- model.estimates$TGR[i] + qt(0.025, model.estimates$DF[i]) * model.estimates$"Std.Error"[i]
            model.estimates$ub95.TGR[i] <- model.estimates$TGR[i] + qt(0.975, model.estimates$DF[i]) * model.estimates$"Std.Error"[i]
            model.estimates$pvalue[i] <- 2 * (1 - pt(abs(model.estimates$"t-value"[i]), model.estimates$DF[i]))
        }
    }
    model.estimates$Contrast <- paste(treatments.numerator, "vs.", treatment.denominator)
    row.names(model.estimates) <- NULL
    # print(model.estimates)
    model.estimates <- model.estimates[, c("Contrast", "TGR", "lb95.TGR", "ub95.TGR", "pvalue")]

    return(model.estimates)
}

tumor.growth.rate.all.therapies.vs.vehicle <- Tumor.growth.curves.exponential(dataset = TGI.dataset, treatments.numerator = c("Therapy_A", "Therapy_B", "Therapy_Combo_AB"), treatment.denominator = "Control")
tumor.growth.rate.therapyA.vs.therapyB <- Tumor.growth.curves.exponential(dataset = TGI.dataset, treatments.numerator = "Therapy_A", treatment.denominator = "Therapy_B", p.tails = 2)

#########################################################################################################################################################################
## Box 6:
## Simulate a data set with quadratic exponential growth
##
## A simulated dataset was created with the following characteristics:
## One vehicle group, 2 single agent therapies (A and B) group, and one combination therapy A+B group ('group.name'). Each group starts with 20 animals ('n.animals').
## Day 0 is when animals are assigned to groups and treatment starts, if any. The study lasts for 24 days ('duration').
## All animals start with a tumor volume of 100 mm^3 ('initial.tumor.volume') and have their tumor volumes measured every 3 days. The measuring error amounts to CV = 25% ('residual.cv').
## The average tumor growth rates (in %) for the groups are: Vehicle = 5; Therapy_A = 4.5; Therapy_B = 3.5; Therapy_AB = 1 ('group.growth.rate.day.perc') with intra-group CV=10% ('growth.rate.cv')
## The average tumor decay rates (in %) for the groups are: Vehicle = 0.2; Therapy_A = 0.125; Therapy_B = 0.1; Therapy_AB = 0 ('group.decay.rate.day.perc') with intra-group CV=10% ('growth.rate.cv')

simulate.quad.exp.growth <- function(group.name, duration, initial.tumor.volume, group.growth.rate.day.perc, growth.rate.cv, group.decay.rate.day.perc, decay.rate.cv, residual.cv, n.animals) {
    group <- list()
    for (i in 1:n.animals) {
        ln.tumor.volume <- numeric()
        ln.tumor.volume[1] <- rnorm(n = 1, mean = log(initial.tumor.volume), sd = sqrt(log((residual.cv^2) + 1)))
        days <- seq(0, duration, 3)
        growth.rate.individual <- rnorm(n = (length(days) - 1), mean = group.growth.rate.day.perc / 100, sd = group.growth.rate.day.perc * growth.rate.cv / 100)
        decay.rate.individual <- rnorm(n = (length(days) - 1), mean = -(group.decay.rate.day.perc / 100), sd = group.decay.rate.day.perc * decay.rate.cv / 100)

        for (t in 2:length(days)) {
            ln.tumor.volume[t] <- rnorm(n = 1, mean = ln.tumor.volume[t - 1] + (growth.rate.individual[t - 1] * days[t]) + (decay.rate.individual[t - 1] * (days[t]^2)), sd = sqrt(log((residual.cv^2) + 1)))
        }
        group[[i]] <- data.frame(Day = days, Animal = paste(group.name, i, sep = "_"), Group = group.name, TV = exp(ln.tumor.volume))
    }
    group <- plyr::rbind.fill(group)
    group
}

set.seed(12345)
vehicle.group.2 <- simulate.quad.exp.growth(
    group.name = "Control", duration = 24, initial.tumor.volume = 100, group.growth.rate.day.perc = 5, growth.rate.cv = 0.1,
    group.decay.rate.day.perc = 0.2, decay.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
therapyA.group.2 <- simulate.quad.exp.growth(
    group.name = "Therapy_A", duration = 24, initial.tumor.volume = 100, group.growth.rate.day.perc = 4.5, growth.rate.cv = 0.1,
    group.decay.rate.day.perc = 0.125, decay.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
therapyB.group.2 <- simulate.quad.exp.growth(
    group.name = "Therapy_B", duration = 24, initial.tumor.volume = 100, group.growth.rate.day.perc = 3.5, growth.rate.cv = 0.1,
    group.decay.rate.day.perc = 0.1, decay.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
therapyAB.group.2 <- simulate.quad.exp.growth(
    group.name = "Therapy_Combo_AB", duration = 24, initial.tumor.volume = 100, group.growth.rate.day.perc = 1, growth.rate.cv = 0.1,
    group.decay.rate.day.perc = 0, decay.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
TGI.dataset.2 <- rbind(vehicle.group.2, therapyA.group.2, therapyB.group.2, therapyAB.group.2)

head(TGI.dataset.2)

TGI.stats.2 <- summary.stats.TGI(dataset = TGI.dataset.2)
plot.TGI.dataset(TGI.stats.2)
plot.TGI.dataset(TGI.stats.2, ln.scale = "yes")

#########################################################################################################################################################################
## Box 7:
## Analysis of growth curves exponential:
##
## Tumor.growth.curves.quadratic.exponential()
## Requires a dataset in the same format as TGI.dataset ('dataset').
## Compares tumor growth rates "TGR" and tumor decay rates "TDR" of treatments ('treatments.numerator') versus control or other treatments ('treatment.denominator').
## One-sided or two-sided tests ('p.tails')

Tumor.growth.curves.quadratic.exponential <- function(dataset, treatments.numerator, treatment.denominator, p.tails = 1) {
    TV.data <- subset(dataset, Group %in% c(treatments.numerator, treatment.denominator))
    TV.data$lnTV <- log(TV.data$TV)
    TV.data$Group <- factor(TV.data$Group, levels = c(treatment.denominator, treatments.numerator))
    TV.data$Day <- as.numeric(TV.data$Day)
    model.fit <- nlme::lme(lnTV ~ Group * poly(Day, 2), random = ~ 1 + poly(Day, 2) | Animal, data = TV.data, na.action = na.omit, method = "REML")
    model.estimates <- as.data.frame(summary(model.fit)$tTable)
	# print(model.estimates)
    model.estimates.1 <- model.estimates[grep(":poly(Day, 2)1", rownames(model.estimates), fixed = TRUE), ]
    model.estimates.2 <- model.estimates[grep(":poly(Day, 2)2", rownames(model.estimates), fixed = TRUE), ]
    names(model.estimates.1)[names(model.estimates.1) == "Value"] <- "TGR"
    model.estimates.1$lb95.TGR <- NA
    model.estimates.1$ub95.TGR <- NA
    model.estimates.1$pval.TGR <- NA
    names(model.estimates.2)[names(model.estimates.2) == "Value"] <- "TDR"
    model.estimates.2$lb95.TDR <- NA
    model.estimates.2$ub95.TDR <- NA
    model.estimates.2$pval.TDR <- NA
    for (i in 1:nrow(model.estimates.1)) {
        if (p.tails == 1) {
            model.estimates.1$ub95.TGR[i] <- model.estimates.1$TGR[i] + qt(0.95, model.estimates.1$DF[i]) * model.estimates.1$"Std.Error"[i]
            model.estimates.1$pval.TGR[i] <- pt(model.estimates.1$"t-value"[i], model.estimates.1$DF[i])
            model.estimates.2$ub95.TDR[i] <- model.estimates.2$TDR[i] + qt(0.95, model.estimates.2$DF[i]) * model.estimates.2$"Std.Error"[i]
            model.estimates.2$pval.TDR[i] <- pt(model.estimates.2$"t-value"[i], model.estimates.2$DF[i])
        }
        if (p.tails == 2) {
            model.estimates.1$lb95.TGR[i] <- model.estimates.1$TGR[i] + qt(0.025, model.estimates.1$DF[i]) * model.estimates.1$"Std.Error"[i]
            model.estimates.1$ub95.TGR[i] <- model.estimates.1$TGR[i] + qt(0.975, model.estimates.1$DF[i]) * model.estimates.1$"Std.Error"[i]
            model.estimates.1$pval.TGR[i] <- 2 * (1 - pt(abs(model.estimates.1$"t-value"[i]), model.estimates.1$DF[i]))
            model.estimates.2$lb95.TDR[i] <- model.estimates.2$TDR[i] + qt(0.025, model.estimates.2$DF[i]) * model.estimates.2$"Std.Error"[i]
            model.estimates.2$ub95.TDR[i] <- model.estimates.2$TDR[i] + qt(0.975, model.estimates.2$DF[i]) * model.estimates.2$"Std.Error"[i]
            model.estimates.2$pval.TDR[i] <- 2 * (1 - pt(abs(model.estimates.2$"t-value"[i]), model.estimates.2$DF[i]))
        }
    }
    model.estimates.1$Contrast <- paste(treatments.numerator, "vs.", treatment.denominator)
    model.estimates.2$Contrast <- paste(treatments.numerator, "vs.", treatment.denominator)
    row.names(model.estimates.1) <- NULL
    row.names(model.estimates.2) <- NULL
    model.estimates <- merge(model.estimates.1[, c("Contrast", "TGR", "lb95.TGR", "ub95.TGR", "pval.TGR")], model.estimates.2[, c("Contrast", "TDR", "lb95.TDR", "ub95.TDR", "pval.TDR")], by = "Contrast", sort = FALSE)
    return(model.estimates)
}

tumor.growth.rate.quad.exp.therapyA.vs.vehicle <- Tumor.growth.curves.quadratic.exponential(dataset = TGI.dataset.2, treatments.numerator = "Therapy_A", treatment.denominator = "Control")
tumor.growth.rate.quad.exp.therapyA.vs.therapyB <- Tumor.growth.curves.quadratic.exponential(dataset = TGI.dataset.2, treatments.numerator = "Therapy_A", treatment.denominator = "Therapy_B", p.tails = 2)


#########################################################################################################################################################################
## Box 8:
## Simulate a data set with gompertzian growth
##
## A simulated dataset was created with the following characteristics:
## One vehicle group, 2 single agent therapies (A and B) group, and one combination therapy A+B group ('group.name'). Each group starts with 20 animals ('n.animals').
## Day 0 is when animals are assigned to groups and treatment starts, if any. The study lasts for 60 days ('duration').
## All animals start with a tumor volume of 100 mm^3 ('initial.tumor.volume') and have their tumor volumes measured every 3 days. The maximum expected a tumor will grow is 6000 mm^3.
## The measuring error amounts to CV = 25% ('residual.cv'). The average tumor growth rates (in %) for the groups are: Vehicle = 6; Therapy_A = 4; Therapy_B = 3; Therapy_AB = 2 ('group.growth.rate.day.perc') 
## with intra-group CV=10% ('growth.rate.cv')

simulate.gompertzian.growth <- function(group.name, duration, initial.tumor.volume, maximum.tumor.volume, group.growth.rate.day.perc, growth.rate.cv, residual.cv, n.animals) {
    group <- list()
    for (i in 1:n.animals) {
        ln.tumor.volume <- numeric()
        ln.tumor.volume[1] <- rnorm(n = 1, mean = log(initial.tumor.volume), sd = sqrt(log((residual.cv^2) + 1)))
        days <- seq(0, duration, 3)
        growth.rate.individual <- rnorm(n = (length(days) - 1), mean = group.growth.rate.day.perc / 100, sd = group.growth.rate.day.perc * growth.rate.cv / 100)
        carrying.capacity <- log(maximum.tumor.volume)

        for (t in 2:length(days)) {
            growth.term <- (ln.tumor.volume[1] - carrying.capacity) * exp(-growth.rate.individual[t - 1] * days[t])
            ln.tumor.volume[t] <- rnorm(n = 1, mean = carrying.capacity + growth.term, sd = sqrt(log((residual.cv^2) + 1)))
        }
        group[[i]] <- data.frame(Day = days, Animal = paste(group.name, i, sep = "_"), Group = group.name, TV = exp(ln.tumor.volume))
    }
    group <- plyr::rbind.fill(group)
    group
}

set.seed(12345)
vehicle.group.3 <- simulate.gompertzian.growth(
    group.name = "Control", duration = 60, initial.tumor.volume = 100, maximum.tumor.volume = 6000,
    group.growth.rate.day.perc = 6, growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
therapyA.group.3 <- simulate.gompertzian.growth(
    group.name = "Therapy_A", duration = 60, initial.tumor.volume = 100, maximum.tumor.volume = 6000,
    group.growth.rate.day.perc = 4, growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
therapyB.group.3 <- simulate.gompertzian.growth(
    group.name = "Therapy_B", duration = 60, initial.tumor.volume = 100, maximum.tumor.volume = 6000,
    group.growth.rate.day.perc = 3, growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
therapyAB.group.3 <- simulate.gompertzian.growth(
    group.name = "Therapy_Combo_AB", duration = 60, initial.tumor.volume = 100, maximum.tumor.volume = 6000,
    group.growth.rate.day.perc = 2, growth.rate.cv = 0.1, residual.cv = 0.25, n.animals = 20
)
TGI.dataset.3 <- rbind(vehicle.group.3, therapyA.group.3, therapyB.group.3, therapyAB.group.3)

head(TGI.dataset.3)

TGI.stats.3 <- summary.stats.TGI(dataset = TGI.dataset.3)
plot.TGI.dataset(TGI.stats.3)
plot.TGI.dataset(TGI.stats.3, ln.scale = "yes")

#########################################################################################################################################################################
## Box 9:
## Analysis of growth curves gompertzian:
##
## Tumor.growth.curves.gompertzian()
## Requires a dataset in the same format as TGI.dataset ('dataset').
## Compares Gompertizian force "Gompertzian.par" of treatments ('treatments.numerator') versus control or other treatments ('treatment.denominator').
## One-sided or two-sided tests ('p.tails')

Tumor.growth.curves.gompertzian <- function(dataset, treatments.numerator, treatment.denominator, p.tails = 1) {
    TV.data <- subset(dataset, Group %in% c(treatments.numerator, treatment.denominator))
    TV.data$y.cross <- log(log(max(TV.data$TV, na.rm = TRUE) / TV.data$TV))
    TV.data$y.cross[!is.finite(TV.data$y.cross)] <- min(TV.data$y.cross[is.finite(TV.data$y.cross)])
    TV.data$Group <- factor(TV.data$Group, levels = c(treatment.denominator, treatments.numerator))
    TV.data$Day <- as.numeric(TV.data$Day)
    model.fit <- nlme::lme(y.cross ~ Group * Day, random = ~ 1 + Day | Animal, data = TV.data, na.action = na.omit, method = "REML", control = nlme::lmeControl(maxIter = 5000, msMaxIter = 5000, niterEM = 2500, msMaxEval = 200))
    model.estimates <- as.data.frame(summary(model.fit)$tTable)
    model.estimates <- model.estimates[grep(":Day", rownames(model.estimates)), ]
    names(model.estimates)[names(model.estimates) == "Value"] <- "Gompertz.par"
    model.estimates$lb95.Gompertz.par <- NA
    model.estimates$ub95.Gompertz.par <- NA
    model.estimates$pval.Gompertz.par <- NA

    for (i in 1:nrow(model.estimates)) {
        if (p.tails == 1) {
            model.estimates$ub95.Gompertz.par[i] <- model.estimates$Gompertz.par[i] + qt(0.95, model.estimates$DF[i]) * model.estimates$"Std.Error"[i]
            model.estimates$pval.Gompertz.par[i] <- 1 - pt(model.estimates$"t-value"[i], model.estimates$DF[i])
        }
        if (p.tails == 2) {
            model.estimates$lb95.Gompertz.par[i] <- model.estimates$Gompertz.par[i] + qt(0.025, model.estimates$DF[i]) * model.estimates$"Std.Error"[i]
            model.estimates$ub95.Gompertz.par[i] <- model.estimates$Gompertz.par[i] + qt(0.975, model.estimates$DF[i]) * model.estimates$"Std.Error"[i]
            model.estimates$pval.Gompertz.par[i] <- 2 * (1 - pt(abs(model.estimates$"t-value"[i]), model.estimates$DF[i]))
        }
    }
    model.estimates$Contrast <- paste(treatments.numerator, "vs.", treatment.denominator)
    row.names(model.estimates) <- NULL
    model.estimates <- model.estimates[, c("Contrast", "Gompertz.par", "lb95.Gompertz.par", "ub95.Gompertz.par", "pval.Gompertz.par")]

    return(model.estimates)
}

tumor.growth.rate.gompertzian.therapyA.vs.vehicle <- Tumor.growth.curves.gompertzian(dataset = TGI.dataset.3, treatments.numerator = c("Therapy_A", "Therapy_B"), treatment.denominator = "Control")
tumor.growth.rate.gompertzian.therapyA.vs.therapyB <- Tumor.growth.curves.gompertzian(dataset = TGI.dataset.3, treatments.numerator = "Therapy_A", treatment.denominator = "Therapy_B", p.tails = 2)


#########################################################################################################################################################################
## Box 10:
## Analysis of area under the curve:
##
## Tumor.AUC()
## Requires a dataset in the same format as TGI.dataset ('dataset').
## Compares the areas under the curves "AUC" from treatments ('treatments.numerator') versus control or other treatments ('treatment.denominator').
## One-sided or two-sided tests ('p.tails')

Tumor.AUC <- function(dataset, treatments.numerator, treatment.denominator, p.tails = 1) {
    TV.data <- subset(dataset, Group %in% c(treatments.numerator, treatment.denominator))
    TV.data$lnTV <- log(TV.data$TV)
    baseline <- min(TV.data$Day)
    TV.baseline <- subset(TV.data, Day %in% baseline)[, c("Animal", "lnTV")]
    names(TV.baseline)[names(TV.baseline) == "lnTV"] <- "lnTV.baseline"
    TV.data <- subset(TV.data, !(Day %in% baseline))
    TV.data <- merge(TV.data, TV.baseline, by = "Animal", sort = FALSE)
    TV.data$height <- TV.data$lnTV - TV.data$lnTV.baseline
    TV.data$height[TV.data$height < 0] <- 0

    AUC.data <- split(TV.data, TV.data$Animal, drop = TRUE)
    AUC.data <- lapply(AUC.data, function(x) {
        new.x <- data.frame(Animal = unique(x$Animal), Group = unique(x$Group))
        x$height <- zoo::na.locf(x$height)
        idx <- order(x$Day)
        new.x$AUC <- sum(diff(x$Day[idx]) * zoo::rollmean(x$height, 2))
        new.x
    })
    AUC.data <- plyr::rbind.fill(AUC.data)
    AUC.data$Group <- factor(AUC.data$Group, levels = c(treatment.denominator, treatments.numerator))

    model.fit <- lm(AUC ~ Group, data = AUC.data)
    model.estimates <- as.data.frame(coef(summary(model.fit)))
    model.estimates <- model.estimates[paste("Group", treatments.numerator, sep = ""), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
    names(model.estimates)[names(model.estimates) == "Estimate"] <- "AUC.diff"

    if (p.tails == 1) {
        model.estimates$lb95.AUC.diff <- NA
        model.estimates$ub95.AUC.diff <- model.estimates$AUC.diff + qt(0.95, model.fit$df.residual) * model.estimates$"Std. Error"
        model.estimates$pvalue <- pt(model.estimates$"t value", model.fit$df.residual)
    }
    if (p.tails == 2) {
        model.estimates$lb95.AUC.diff <- model.estimates$AUC.diff + qt(0.025, model.fit$df.residual) * model.estimates$"Std. Error"
        model.estimates$ub95.AUC.diff <- model.estimates$AUC.diff + qt(0.975, model.fit$df.residual) * model.estimates$"Std. Error"
        model.estimates$pvalue <- 2 * (1 - pt(abs(model.estimates$"t value"), model.fit$df.residual))
    }
    model.estimates$Contrast <- paste(treatments.numerator, "vs.", treatment.denominator)
    row.names(model.estimates) <- NULL
    model.estimates <- model.estimates[, c("Contrast", "AUC.diff", "lb95.AUC.diff", "ub95.AUC.diff", "pvalue")]

    return(model.estimates)
}

AUC.all.therapies.vs.vehicle <- Tumor.AUC(dataset = TGI.dataset, treatments.numerator = c("Therapy_A", "Therapy_B", "Therapy_Combo_AB"), treatment.denominator = "Control")
AUC.therapyA.vs.therapyB <- Tumor.AUC(dataset = TGI.dataset, treatments.numerator = "Therapy_A", treatment.denominator = "Therapy_B", p.tails = 2)

#########################################################################################################################################################################
## Box 11:
## TTE:
##
## TTE.analysis()
## Requires a dataset in the same format as TGI.dataset ('dataset'). Requires an event based on tumor volumes values ('event')
## Compares treatments ('treatments.numerator') to control ('treatment.denominator').
## Calculates results for Cox-Proportional Hazards regression ('method' = "cox") or log-rank tests ('method' = "log-rank")
## If 'plotKM' = TRUE, Kaplan-Meier curves are plotted in colors (default) or in gray-shades (if grey.shades = "yes")

TTE.analysis <- function(dataset, event, treatments.numerator, treatment.denominator, method, plotKM = FALSE, grey.shades = "no"){
	dataset         <- subset(dataset, Group %in% c(treatments.numerator, treatment.denominator))
	dataset$Day     <- as.numeric(dataset$Day)
	day0            <- min(dataset$Day)
	day0.tv         <- dataset[dataset$Day == day0, "TV"]
	data.n          <- tapply(dataset$TV, dataset$Animal, length)
	dataset$day0.TV <- rep(day0.tv, data.n)
	dataset         <- dataset[!is.na(dataset$TV), ]
	dataset$event   <- dataset$TV > event
	
	data.tte   <- list()
	animal.ids <- unique(dataset$Animal)
	for(i in animal.ids){
		surv.time1 <- surv.time2 <- NA
		data.tte1  <- dataset[dataset$Animal == i, ]
		if(all(!data.tte1$event, na.rm = T)){
			surv.time2 <- max(data.tte1[!is.na(data.tte1$TV), "Day"], na.rm = T) 
			event.day  <- which(data.tte1$Day == surv.time2)
		} else {
			surv.time1 <- min(data.tte1[data.tte1$event, "Day"], na.rm = T)
			event.day  <- which(data.tte1$Day == surv.time1) 
		}
		data.tte <- rbind(data.tte, cbind(data.tte1[event.day, ], surv.time = max(c(surv.time1,surv.time2), na.rm = T) - day0))
	}
  
	data.tte$Group <- factor(data.tte$Group, levels = c(treatment.denominator, treatments.numerator))

	if(method == "cox"){
		cox              <- survival::coxph(survival::Surv(surv.time,event) ~ Group, data = data.tte)
		res              <- cbind(summary(cox)$coef, summary(cox)$conf.int[, 3:4, drop = F])
		rownames(res)    <- gsub("Group", "", rownames(res))
		colnames(res)[colnames(res) == "Pr(>|z|)"] <- "pvalue"
		colnames(res)    <- gsub(".95", "HR", colnames(res))
		colnames(res)[2] <- "HR"
		res              <- data.frame(Contrast = paste(rownames(res), "vs.", cox$xlevels$Group[1]), res)
		rownames(res)    <- NULL
	}

	if(method == "log-rank"){
		res <- list()
		for(j in 1:length(treatments.numerator)){
			lr  <- survival::survdiff(survival::Surv(surv.time, event) ~ Group, data = data.tte, Group %in% c(treatments.numerator[j], treatment.denominator))
			res <- rbind(res, data.frame(Contrast = paste(treatments.numerator[j], "vs.", treatment.denominator), chisq = lr$chisq, pvalue = 1 - pchisq(lr$chisq, 1)))
		}
	}
	
	if(plotKM){
		surv.fit <- survival::survfit(survival::Surv(surv.time, event) ~ Group, data = data.tte)
		grps     <- gsub("Group=", "", names(surv.fit$strata))
		grps.col <- viridis::viridis(length(grps))
		if(grey.shades == "yes") grps.col <- grDevices::gray.colors(n = length(grps))
		plot(surv.fit, conf.int = FALSE, mark.time = TRUE, lwd = 2, ylab = " 1 - Pr(event)", xlab = "Days", main = "Kaplan-Meier curves", cex.lab = 1.3, font.lab = 2, cex.main = 2.5, font.main = 2, col = grps.col, xaxt = "n", yaxt = "n")
		abline(v = unique(dataset$Day), col = "grey90", lty = 3, lwd = 0.7)
		abline(h = c(0, 0.25, 0.5, 0.75, 1), col = "grey90", lty = 3, lwd = 0.7)
		axis(1, at = unique(dataset$Day), label = unique(dataset$Day))
		axis(2, at = c(0, 0.25, 0.5, 0.75, 1), label = c(0, 0.25, 0.5, 0.75, 1), las = 2)
		legend("bottomleft", legend = grps, col = grps.col, bty = "n", lwd = 1, cex = 1.2)
		lines(surv.fit, conf.int = FALSE, mark.time = TRUE, col = grps.col)
	}
	return(res)
}

tte.cox.therapies.A.and.B.vs.control <- TTE.analysis(TGI.dataset.3, event = 2000, treatments.numerator = c("Therapy_A", "Therapy_B"), treatment.denominator = "Control", method = "cox", plotKM = FALSE)
tte.lr.therapies.A.and.B.vs.control  <- TTE.analysis(TGI.dataset.3, event = 2000, treatments.numerator = c("Therapy_A", "Therapy_B"), treatment.denominator = "Control", method = "log-rank", plotKM = TRUE)

#########################################################################################################################################################################
## Box 12:
## Categorical:
##
## Response.analysis()
## Requires a dataset in the same format as TGI.dataset ('dataset'). Compares treatments ('treatments.numerator') to control ('treatment.denominator').
## Requires the initial (day.0) and final ('day.resp') days or comparison
## PD threshold ('PD.cutoff') above which is PD, that is, PD = T/C > 1 + PD.cutoff 
## PR threshold ('PR.cutoff') below which is PR, that is, PR = T/C < 1 - PR.cutoff
## Should SD be included ('SD.include') as nonresponder (1) or responder (2) or excluded (0)
## One-sided or two-sided tests ('p.tails')

Response.analysis <- function(dataset, treatments.numerator, treatment.denominator, day.0, day.resp, PD.cutoff = 0.5, PR.cutoff = 0.3, SD.include = 1, p.tails = 2){
	TV.data.day0  <- subset(dataset, Day == day.0 & Group %in% c(treatments.numerator, treatment.denominator))
	TV.data.final <- subset(dataset, Day == day.resp & Group %in% c(treatments.numerator, treatment.denominator))
	response.data <- merge(TV.data.day0, TV.data.final[, c("Animal", "TV")], by = "Animal", suffix = c("", ".final"), sort = FALSE)
	
	response.data$resp4 <- (response.data[, "TV.final"] / response.data[, "TV"] >= (1 + PD.cutoff)) + 
						   (response.data[, "TV.final"] / response.data[, "TV"] < (1 + PD.cutoff) & response.data[, "TV.final"] / response.data[, "TV"] > (1 - PR.cutoff)) * 2 + 
						   (response.data[, "TV.final"] / response.data[, "TV"] < (1 - PR.cutoff) & response.data[, "TV.final"] > 10) * 3 +
						   (response.data[, "TV.final"] <= 10) * 4
    response.data$resp4[is.na(response.data[, "TV.final"])] <- NA
	response.data$resp=NA
	if (SD.include > 0){
		response.data$resp <-(response.data$resp4 >= (4 - SD.include)) + 0 
	} else {
		response.data$resp[response.data$resp4 == 1] = 0
		response.data$resp[response.data$resp4 > 2] = 1
	}
	response.data$resp4 <- factor(response.data$resp4, levels = 1:4)
	response.data$resp  <- factor(response.data$resp, levels = c(0,1))
	rt4                 <- table(response.data$Group, response.data$resp4)
	rt2                 <- table(response.data$Group, response.data$resp)
	colnames(rt4)       <- c("PD", "SD", "PR", "CR")
	colnames(rt2)       <- c("NonReponder", "Responder")
 
	resp.pv <- p.fisher <- rep(1, nrow(rt4) - 1)
	for(i in 2:nrow(rt4)){ 
		if(sum(apply(rt4[c(1, i), ], 2, sum) == 0) < 3){
			resp.pv[i - 1] <- coin::pvalue(coin::independence_test(resp4 ~ as.factor(Group), data = response.data[response.data$Group %in% rownames(rt4)[c(1, i)], ],
																   distribution = "exact", ytrafo = function(x) coin::rank_trafo(x, ties.method = c("mid-ranks", "random"))))/((p.tails == 1) + 1)
			p.fisher[i-1]  <- fisher.test(rt2[c(1,i),])$p.value/((p.tails == 1) + 1) 
		}
	}
	
	res <- data.frame(Contrast = paste(rownames(rt4)[-1], "vs.", rownames(rt4)[1]), p.trend = resp.pv, p.fisher = p.fisher)
 
	return(res)
}

#########################################################################################################################################################################
## Box 13:
## Function to plot tumor growth curves with Bliss line:
##
## plot.Bliss()
## Requires a dataset in the same format as TGI.dataset ('dataset'). 
## The names of control ('control.name'), single treatments ('single.1.name' and 'single.2.name'), and combo ('combo.name') groups are required.
## Depicts a plot with tumor growth curves in the original scale (default) or in natural ln-scale (if ln.scale = "yes")
## Tumor growth curves are plotted in colors (default) or in gray-shades (if grey.shades = "yes")

plot.Bliss <- function(dataset, control.name, single.1.name, single.2.name, combo.name, ln.scale = "no", grey.shades = "no"){
	dataset.list <- subset(dataset, Group %in% c(control.name, single.1.name, single.2.name, combo.name))
	dataset.list <- split(dataset.list, list(dataset.list$Group, dataset.list$Day, drop = TRUE))
	dataset.list <- lapply(dataset.list, function(x){
						new.x          <- as.data.frame(unique(x[, c("Day", "Group")]))
						new.x$est      <- exp(mean(log(x$TV), na.rm = TRUE))
						new.x$se       <- sqrt(var(log(x$TV), na.rm = TRUE)/nrow(x[!is.na(x$TV), ]))
						new.x$lb95.est <- exp(log(new.x$est) + qnorm(0.025) * new.x$se)
						new.x$ub95.est <- exp(log(new.x$est) + qnorm(0.975) * new.x$se)
						new.x
					})
	dataset.list <- plyr::rbind.fill(dataset.list)
	
	dataset.list <- split(dataset.list, list(dataset.list$Day), drop = TRUE)
	dataset.list <- lapply(dataset.list, function(x){
						new.x     <- data.frame(Day = unique(x$Day), Group = "Bliss")
						new.x$est <- x$est[x$Group == single.1.name] * x$est[x$Group == single.2.name] / x$est[x$Group == control.name]
						plyr::rbind.fill(x, new.x)
					})
	summ.stats   <- plyr::rbind.fill(dataset.list)

	summ.stats$Group       <- factor(summ.stats$Group, levels = c(control.name, single.1.name, single.2.name, combo.name, "Bliss"))
	summ.stats$group.color <- viridis::viridis(nlevels(summ.stats$Group) + 1)[as.numeric(summ.stats$Group)] ## requires package 'viridis' to be installed
	summ.stats$group.color[summ.stats$Group == "Bliss"] <- "#FF0000FF"
	if(grey.shades == "yes") {
		summ.stats$group.color <- grDevices::gray.colors(n = nlevels(summ.stats$Group) + 1)[as.numeric(summ.stats$Group)] ## requires package 'grDevices' to be installed
		summ.stats$group.color[summ.stats$Group == "Bliss"] <- "#000000FF"
	}
	summ.stats$group.pch   <- c(15, 16, 17, 18, NA)[as.numeric(summ.stats$Group)]
	y.ticks                <- seq(0, max(summ.stats$ub95.est, na.rm = TRUE), 500)
	
	if(ln.scale == "yes"){
		summ.stats$est      <- log(summ.stats$est)
		summ.stats$lb95.est <- log(summ.stats$lb95.est)
		summ.stats$ub95.est <- log(summ.stats$ub95.est)
		y.ticks             <- log(y.ticks)
	}
	par(mar = c(5, 6, 3, 1))
	plotrix::plotCI(x = summ.stats$Day, y = summ.stats$est, li = summ.stats$lb95.est, ui = summ.stats$ub95.est, xaxt = "n", xlab = "Days", yaxt = "n", ylab = "", cex.lab = 1.3, font.lab = 2, main = "Tumor Growth Curves", cex.main = 2.5, font.main = 2, pch = NA, slty = 0, lwd = 1.7)
	axis(1, at = unique(summ.stats$Day), label = unique(summ.stats$Day))
	abline(v = unique(summ.stats$Day), lty = 3, col = "gray80")
	if(ln.scale == "yes"){
		axis(2, at = y.ticks, label = format(exp(y.ticks), sci = T), las = 2)
	} else {
		axis(2, at = y.ticks, label = format(y.ticks, sci = T), las = 2)
	}
	abline(h = y.ticks, lty = 3, col = "gray80")
	mtext(side = 2, text = "Tumor Volume (mm^3)", line = 5, font = 2, cex = 1.3)
	group.data <- subset(summ.stats, Group != "Bliss")
	group.data <- split(group.data, group.data$Group, drop = TRUE)
	invisible(lapply(group.data, function(x){
					plotrix::plotCI(x = x$Day, y = x$est, li = x$lb95.est, ui = x$ub95.est, pch = unique(x$group.pch), col = unique(x$group.color), add = TRUE, cex = 1.5)
					lines(x = x$Day, y = x$est, col = unique(x$group.color), lwd = 2)
			 }))
	bliss.line <- subset(summ.stats, Group == "Bliss")
	lines(x = bliss.line$Day, y = bliss.line$est, col = unique(bliss.line$group.color), lwd = 2, lty = 3)
	
	legend.attributes <- unique(summ.stats[, c("Group", "group.color", "group.pch")])
	legend("topleft", legend = legend.attributes$Group, col = legend.attributes$group.color, pch = legend.attributes$group.pch, lwd = 1.4, lty = c(1, 1, 1, 1, 3), cex = 1.3, bty = "n")
}

plot.Bliss(TGI.dataset, control.name = "Control", single.1.name = "Therapy_A", single.2.name = "Therapy_B", combo.name = "Therapy_Combo_AB")
plot.Bliss(TGI.dataset, control.name = "Control", single.1.name = "Therapy_A", single.2.name = "Therapy_B", combo.name = "Therapy_Combo_AB", ln.scale = "yes")
plot.Bliss(TGI.dataset, control.name = "Control", single.1.name = "Therapy_A", single.2.name = "Therapy_B", combo.name = "Therapy_Combo_AB", ln.scale = "yes", grey.shades = "yes")

#########################################################################################################################################################################
## Box 14:
## Bliss independence test for synergy:
##
## Bliss.synergy.test()
## Requires a dataset in the same format as TGI.dataset ('dataset'). Estimates TGI for multiple time-points ('time.points').
## TGI estimates can be either adjusted or unadjusted ('baseline.adjusted') for baseline day ('baseline'). 
## The names of control ('control.name'), single treatments ('single.1.name' and 'single.2.name'), and combo ('combo.name') groups are required.

Bliss.synergy.test <- function(dataset, time.points, baseline.adjusted = "no", baseline = NULL, control.name, single.1.name, single.2.name, combo.name){

	if(baseline.adjusted == "yes" & is.null(baseline)) stop("Please either assign a baseline time point for TGI calculations or set 'baseline.adjusted' to 'no'.")
	TV.data      <- subset(dataset, Day %in% c(baseline, time.points) & Group %in% c(control.name, single.1.name, single.2.name, combo.name))
	TV.data$lnTV <- log(TV.data$TV)
	if(baseline.adjusted == "yes"){
			TV.baseline <- subset(TV.data, Day %in% baseline)[, c("Animal", "lnTV")]
			names(TV.baseline)[names(TV.baseline) == "lnTV"] <- "lnTV.baseline"
			TV.data     <- subset(TV.data, !(Day %in% baseline))
			TV.data     <- merge(TV.data, TV.baseline, by = "Animal", sort = FALSE)
	}	
	TV.data$Group <- factor(TV.data$Group, levels = c(control.name, single.1.name, single.2.name, combo.name))
	day.idx       <- sort(unique(TV.data$Day))
	TV.data$Day   <- factor(TV.data$Day, levels = day.idx)

	model.estimates <- list()
	for(i in 1:length(day.idx)){
		TV.data.day          <- subset(TV.data, Day == day.idx[i])
		if(baseline.adjusted == "no")  model.fit <- lm(lnTV ~ Group, data = TV.data.day)
		if(baseline.adjusted == "YES") model.fit <- lm(lnTV ~ lnTV.baseline + Group, data = TV.data.day)
		fit.combo            <- gmodels::fit.contrast(model.fit, "Group", c(1, -1, -1, 1), conf.int = 0.9, df = TRUE)
		model.estimates[[i]] <- data.frame(Contrast     = paste(combo.name, " vs. Bliss", sep = ""),
										   Day          = day.idx[i],
										   Synergy      = exp(as.numeric(as.vector(fit.combo[, "Estimate"]))),
										   lb95.Synergy = exp(as.numeric(as.vector(fit.combo[, "lower CI"]))),
										   ub95.Synergy = exp(as.numeric(as.vector(fit.combo[, "upper CI"]))),
										   pvalue       = pt(fit.combo[, "t value"], df = fit.combo[, "DF"]))
	}
	model.estimates <- plyr::rbind.fill(model.estimates)

	return(model.estimates)
}

Bliss.synergy.test.output.unadjusted <- Bliss.synergy.test(dataset = TGI.dataset, time.points = c(0, 3, 6, 9, 12, 15, 18), baseline.adjusted = "no", baseline = NULL, 
														   control.name = "Control", single.1.name = "Therapy_A", single.2.name = "Therapy_B", combo.name = "Therapy_Combo_AB")
Bliss.synergy.test.output.adjusted   <- Bliss.synergy.test(dataset = TGI.dataset, time.points = c(0, 3, 6, 9, 12, 15, 18), baseline.adjusted = "yes", baseline = 0, 
														   control.name = "Control", single.1.name = "Therapy_A", single.2.name = "Therapy_B", combo.name = "Therapy_Combo_AB")


#########################################################################################################################################################################
## Box 15:
## Sample size estimation:
##
## Power.function()
## Number of groups ('ngrp'), one-sided or two-sided ('side'), or significance level ('alpha').
## Either 'method' = "ANOVA" if no additional covariates are expected, or 'method' = "ANCOVA" when covariates are present.
## One and only one of these need to be NULL: residual standard deviation ('sigma'); sample size ('n'); effect size ('eff'); power ('power').

pwr.gt2.grps <- function(n, sigma, eff, side, ngrp, alpha, method = "ANOVA"){ 
	DFn    <- (n - 1) * ngrp - (method == "ANCOVA")
	SEn    <- sqrt(2/n) * sigma
	tcrit  <- qt(1 - alpha/side, DFn)
	noncen <- log(1 + eff)/SEn
	1 - pt(tcrit, DFn, ncp = noncen) + pt(-tcrit, DFn, ncp = noncen)
}

pwr.2.grps <- function(n1, n2, sigma1, sigma2, eff, side, alpha, method = "ANOVA"){ 
	DFn    <- n1 + n2 - 2 - (method == "ANCOVA")
    SEn    <- sqrt(1/n1 * sigma1^2 + 1/n2 * sigma2^2)
	tcrit  <- qt(1 - alpha/side, DFn)
	noncen <- log(1 + eff)/SEn
	1 - pt(tcrit, DFn, ncp = noncen) + pt(-tcrit, DFn, ncp = noncen)
}

Power.function <- function(sigma = 1, n = NULL, eff = .5, power = 0.8, ngrp = 2, alpha = 0.05, side = 2, method = "ANCOVA"){ 
	pars  <- c(is.null(sigma),is.null(n),is.null(eff),is.null(power))
    if (sum(pars) != 1) return("Error: one and only one parameter should be NULL") 
    par   <- which(pars)
    res.p <- list()

    if (par == 4 & ngrp > 2){
		for (ns in n){
			for (s in sigma){
				for (ef in eff){
					DFn   <- (ns - 1) * ngrp - (method == "ANCOVA")
					SEn   <- sqrt(2/ns) * s
					pow   <- pwr.gt2.grps(n = ns, sigma = s, eff = ef, side = side, ngrp = ngrp, alpha = alpha, method = method)
					res1  <- data.frame(sigma = s, ngroup = ngrp, N = ns, DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
					res.p <- rbind(res.p, res1)
				}
			}
		}
	}else{
		if(par == 4 & ngrp == 2){
			if(is.list(sigma)){
				s1 <- sigma[[1]]
				s2 <- sigma[[2]]
			}else{
				s1 <- s2 <- sigma
			}
			if(is.list(n)){
				n1 <- n[[1]]
				n2 <- n[[2]]
			}else{
				n1 <- n2 <- n
			}
			for(ni in 1:length(n1)){
				for(si in 1:length(s1)){
					for (ef in eff){
						DFn   <- (n1[ni] + n2[ni] - 2) - (method == "ANCOVA")
						SEn   <- sqrt(1/n1[ni] * s1[si]^2 + 1/n2[ni] * s2[si]^2)
						pow   <- pwr.2.grps(n1 = n1[ni], n2 = n2[ni], sigma1 = s1[si], sigma2 = s2[si], eff = ef, side = side, alpha = alpha, method = method)
						res1  <- data.frame(sigma1 = s1[si], sigma2 = s2[si], n1 = n1[ni], n2 = n2[ni], DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
						res.p <- rbind(res.p,res1)
					}
				}
			}
		}
	}
   
	if(par != 4 & ngrp > 2){
		for(pow in power){
			if(is.null(sigma)){
				for(ns in n){
					for(ef in eff){
						f     <- function(s) pow - pwr.gt2.grps(n = ns, sigma = s, eff = ef, side = side, ngrp = ngrp, alpha = alpha, method = method)
						sc    <- uniroot(f, c(0.01, 2))$root
						DFn   <- (ns - 1) * ngrp - (method == "ANCOVA")
						SEn   <- sqrt(2/ns) * sc
						res1  <- data.frame(sigma = sc, ngroup = ngrp, N = ns, DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
						res.p <- rbind(res.p, res1)
					}
				}
			}else{
				if(is.null(n)){ 
					for(s in sigma){
						for (ef in eff){
							f     <- function(n) pow - pwr.gt2.grps(n = n, sigma = s, eff = ef, side = side, ngrp = ngrp, alpha = alpha, method = method)
							ns    <- uniroot(f, c(2, 2000))$root
							DFn   <- (ns - 1) * ngrp - (method == "ANCOVA")
							SEn   <- sqrt(2/ns) * s
							res1  <- data.frame(sigma = s, ngroup = ngrp, N = ns, DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
							res.p <- rbind(res.p, res1)
						}
					}
				}else{
					if(is.null(eff)){
						for(s in sigma){
							for(ns in n){
								f     <- function(ef) pow - pwr.gt2.grps(n = ns, sigma = s, eff = ef, side = side, ngrp = ngrp, alpha = alpha, method = method)
								ef    <- uniroot(f, c(0.1,10))$root
								DFn   <- (ns - 1) * ngrp - (method == "ANCOVA")
								SEn   <- sqrt(2/ns) * s
								res1  <- data.frame(sigma = s, ngroup = ngrp, N = ns, DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
								res.p <- rbind(res.p, res1)
							}
						}
					}
				}
			}
		}
	}

	if(par != 4 & ngrp == 2){
		for(pow in power){
			if(is.null(sigma)){
				if(is.list(n)){
					n1 <- n[[1]]
					n2 <- n[[2]]
				} else {
					n1 <- n2 <- n
				}
				
				for(ni in 1:length(n1)){
					for(ef in eff){
						f     <- function(s) pow - pwr.2.grps(n1 = n1[ni], n2 = n2[ni], sigma1 = s, sigma2 = s, eff = ef, side = side, ngrp = ngrp, alpha = alpha, method = method)
						sc    <- uniroot(f, c(0.01, 2))$root
						DFn   <- (n1[ni] + n2[ni] - 2) - (method == "ANCOVA")
						SEn   <- sqrt(1/n1[ni] * sc^2 + 1/n2[ni] * sc^2)
						res1  <- data.frame(sigma1 = sc, sigma2 = sc, n1 = n1[ni], n2 = n2[ni], DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
						res.p <- rbind(res.p, res1)
					}
				}
			}else{
				if(is.null(n)){
					if(is.list(sigma)){
						s1 <- sigma[[1]]
						s2 <- sigma[[2]]
					}else{
						s1 <- s2 <- sigma
					}
					for(si in 1:length(s1)){
						for(ef in eff){
							f     <- function(n) pow - pwr.2.grps(n1 = n, n2 = n, sigma1 = s1[si], sigma2 = s2[si], eff = ef, side = side, alpha = alpha, method = method)
							ns    <- uniroot(f, c(2, 2000))$root
							DFn   <- (ns + ns - 2) - (method == "ANCOVA")
							SEn   <- sqrt(1/ns * s1[si]^2 + 1/ns * s2[si]^2)
							res1  <- data.frame(sigma1 = s1[si], sigma2 = s2[si], n1 = ns, n2 = ns, DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
							res.p <- rbind(res.p, res1)
						}
					}
				}else{
					if(is.null(eff)){
						if(is.list(sigma)){
							s1 <- sigma[[1]]
							s2 <- sigma[[2]]
						}else{
							s1 <- s2 <- sigma
						}
						if(is.list(n)){
							n1 <- n[[1]]
							n2 <- n[[2]]
						}else{
							n1 <- n2 <- n
						}
						for(si in 1:length(s1)){
							for(ni in 1:length(n1)){
								f     <- function(ef) pow - pwr.2.grps(n1 = n1[ni], n2 = n2[ni], sigma1 = s1[si], sigma2 = s2[si], eff = ef, side = side, alpha = alpha, method = method)
								ef    <- uniroot(f, c(0.1, 10))$root
								DFn   <- (n1[ni] + n2[ni] - 2) - (method == "ANCOVA")
								SEn   <- sqrt(1/n1[ni] * s1[si]^2 + 1/n2[ni] * s2[si]^2)
								res1  <- data.frame(sigma1 = s1[si], sigma2 = s2[si], n1 = n1[ni], n2 = n2[ni], DF = DFn, StdErr = SEn, effectsize = ef, power = pow)
								res.p <- rbind(res.p, res1)
							}
						}
					}
				}
			}
		}
	}

	res.p
}

sample.size.ancova <- Power.function(sigma = 0.25, n = NULL, eff = 1, power = 0.8, ngrp = 4, alpha = 0.05, side = 2, method = "ANCOVA")
sample.size.anova  <- Power.function(sigma = 0.25, n = NULL, eff = 1, power = 0.8, ngrp = 4, alpha = 0.05, side = 2, method = "ANOVA")



