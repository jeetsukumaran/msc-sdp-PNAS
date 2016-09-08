require(ggplot2)
require(grid)
require(plyr)

statistical.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

analyze.num.orthospecies.table = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    ddply(d, ~speciation_completion_rate_level, summarise,
        max.lineages=max(num_extant_lineages),
        min.lineages=min(num_extant_lineages),
        median.lineages=median(num_extant_lineages),
        mean.lineages=mean(num_extant_lineages),
        statistical.mode.lineages=statistical.mode(num_extant_lineages),
        sd.lineages=sd(num_extant_lineages),
        max.true.species=max(num_extant_orthospecies),
        min.true.species=min(num_extant_orthospecies),
        median.orthospecies=median(num_extant_orthospecies),
        mean.orthospecies=mean(num_extant_orthospecies),
        statistical.mode.orthospecies=statistical.mode(num_extant_orthospecies),
        sd.orthospecies=sd(num_extant_orthospecies)
        )
}

analyze.error.table.for.inferred.species.vs.true.species = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    d$absolute.error = (d$num_estimated_orthospecies - d$num_extant_orthospecies)
    d$normalized.error = (d$num_estimated_orthospecies - d$num_extant_orthospecies) / d$num_extant_orthospecies
    d$rmse = sqrt(  (d$num_estimated_orthospecies - d$num_extant_orthospecies)  ^ 2 )
    ddply(d, ~speciation_completion_rate_level, summarise,
        max.absolute.error=max(absolute.error),
        min.absolute.error=min(absolute.error),
        mean.absolute.error=mean(absolute.error),
        mean.normalized.error=mean(normalized.error),
        mean.rmse=mean(rmse)
        )
}

analyze.error.table.for.inferred.species.vs.lineages = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    d$absolute.error = (d$num_estimated_orthospecies - d$num_extant_lineages)
    d$normalized.error = (d$num_estimated_orthospecies - d$num_extant_lineages) / d$num_extant_lineages
    d$rmse = sqrt(  (d$num_estimated_orthospecies - d$num_extant_lineages)  ^ 2 )
    ddply(d, ~speciation_completion_rate_level, summarise,
        max.absolute.error=max(absolute.error),
        min.absolute.error=min(absolute.error),
        mean.absolute.error=mean(absolute.error),
        mean.normalized.error=mean(normalized.error),
        mean.rmse=mean(rmse)
        )
}

plot.num.orthospecies.by.speciation.conversion.rate = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    p = ggplot(d, aes(num_extant_orthospecies))
    p = p + geom_histogram(binwidth=1)
    p = p + facet_grid(~speciation_completion_rate_level)
    p
}

plot.inferred.species.vs.true.species = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    maxspp = max(d$num_extant_orthospecies, d$num_estimated_orthospecies)
    p = ggplot(d, aes(num_extant_orthospecies, num_estimated_orthospecies))
    p = p + geom_jitter(aes(color=speciation_completion_rate_level))
    p = p + scale_colour_discrete(name="Conversion Rate")
    p = p + xlab("Simulated Number of Species")
    p = p + ylab("Inferred Number of Species")
    p = p + scale_x_continuous(limits=c(0, maxspp))
    p = p + scale_y_continuous(limits=c(0, maxspp))
    p = p + theme(legend.position = "bottom")
    p = p + geom_abline(intercept = 0, slope = 1, linetype="dotted")
    p = p + theme(
            legend.key.size=unit(1, "cm"),
            legend.title=element_text(size=18),
            legend.text=element_text(size=18),
            axis.text=element_text(size=16,
            face="bold"),
            axis.title=element_text(size=18,
            face="bold"))
    pdf("inferred-ntax-vs-true-ntax.pdf")
    print(p)
    dev.off()
    p
}

plot.inferred.species.vs.lineages = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    maxspp = max(d$num_extant_orthospecies, d$num_estimated_orthospecies)
    p = ggplot(d, aes(num_extant_lineages, num_estimated_orthospecies))
    p = p + geom_jitter(aes(color=speciation_completion_rate_level))
    p = p + scale_colour_discrete(name="Conversion Rate")
    p = p + xlab("Simulated Number of Lineages")
    p = p + ylab("Inferred Number of Species")
    p = p + scale_x_continuous(limits=c(0, maxspp))
    p = p + scale_y_continuous(limits=c(0, maxspp))
    p = p + theme(legend.position = "bottom")
    p = p + geom_abline(intercept = 0, slope = 1, linetype="dotted")
    p = p + theme(
            legend.key.size=unit(1, "cm"),
            legend.title=element_text(size=18),
            legend.text=element_text(size=18),
            axis.text=element_text(size=16,
            face="bold"),
            axis.title=element_text(size=18,
            face="bold"))
    pdf("inferred-ntax-vs-nlineages.pdf")
    print(p)
    dev.off()
    p
}

analyze.error = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    d$error = (d$num_estimated_orthospecies - d$num_extant_orthospecies) / d$num_extant_orthospecies
    ggplot(d, aes(speciation_completion_rate_level, error)) + geom_point()
}

analyze.rmse = function(d) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    d$rmse = sqrt(  (d$num_estimated_orthospecies - d$num_extant_orthospecies)  ^ 2 )
    ggplot(d, aes(speciation_completion_rate_level, rmse)) + geom_point()
}

analyze.fixed.number.of.species.orthospecies.histograms = function(d, true.num.species) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    p = ggplot(d, aes(num_estimated_orthospecies, fill=speciation_completion_rate_level)) + geom_density(alpha=0.2)
    p = p + scale_x_continuous(limits=c(0, 20))
    p = p + scale_y_continuous(limits=c(0, 0.34))
    p = p + geom_vline(aes_string(xintercept=5), color="blue", linetype="dotted", size=0.75)
    p = p + geom_text(data = NULL, x = 9, y = 0.34, label = " <- True number of species = 8", size=5, color="blue")
    p = p + theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold"))
    print(p)
    p
}

plot.bpp.performance.boxplot = function(d, true.num.species) {
    d$speciation_completion_rate_level = factor(d$speciation_completion_rate)
    p = ggplot(d, aes(x=speciation_completion_rate_level, y=num_estimated_orthospecies))
    p = p + geom_boxplot()
    # p = p + geom_jitter(position = position_jitter(height=.1, width=.1))
    p = p + scale_y_continuous(limits=c(0, 20))

    p = p + theme(
            # axis.text.x=element_blank(),
            # axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()#,
            # axis.line=element_blank(),
            # axis.ticks=element_blank(),
            # legend.position="none",
            # panel.background=element_blank(),
            # panel.border=element_blank(),
            # panel.grid.major=element_blank(),
            # panel.grid.minor=element_blank(),
            # plot.background=element_blank()
            )

    p = p + theme(axis.text=element_text(size=16,face="bold"), axis.title=element_text(size=14,face="bold"))
    p = p + geom_hline(aes_string(yintercept=true.num.species), color="blue", linetype="dashed", size=0.75)
    # p = p + ylab("Inferred Number of Species")
    # p = p + xlab("Conversion Rate")
    # p = p + geom_text(data = NULL, x = 1.5, y = true.num.species - .5, label = paste("True number of species =", true.num.species), size=5, color="blue")
    # p = p + theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold"))

    pdf(paste("fixed-number-of-species-performance-spp", true.num.species, ".pdf", sep=""))
    print(p)
    dev.off()
    p
}

# plot.bpp.performance.violin = function(d, true.num.species) {
#     d$conversion.rate.level = factor(d$conversion.rate)
#     p = ggplot(d, aes(x=conversion.rate.level, y=estimated.number.of.species)) + geom_violin(trim=T)
#     p = p + scale_y_continuous(limits=c(0, 20))
#     p = p + geom_hline(aes_string(yintercept=5), color="blue", linetype="dashed", size=0.75)
#     p = p + geom_text(data = NULL, x = 1.5, y = 4.5, label = "True number of species = 5", size=5, color="blue")
#     p = p + theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold"))
# }
# plot.bpp.performance.scatter = function(d, true.num.species) {
#     d$conversion.rate.level = factor(d$conversion.rate)
#     p = ggplot(d)
#     p = p + geom_point(aes(x=conversion.rate.level, y=estimated.number.of.species))
#     p = p + scale_y_continuous(limits=c(0, 20))
#     p = p + geom_hline(aes_string(yintercept=5), color="blue", linetype="dashed", size=0.75)
#     p = p + geom_text(data = NULL, x = 1.5, y = 4.5, label = "True number of species = 5", size=5, color="blue")
#     p = p + theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold"))
# }
# plot.bpp.performance.jitter = function(d, true.num.species) {
#     d$conversion.rate.level = factor(d$conversion.rate)
#     p = ggplot(d)
#     p = p + geom_jitter(aes(x=conversion.rate.level, y=estimated.number.of.species), position = position_jitter(height=.1, width=.1))
#     p = p + scale_y_continuous(limits=c(0, 20))
#     p = p + geom_hline(aes_string(yintercept=5), color="blue", linetype="dashed", size=0.75)
#     p = p + geom_text(data = NULL, x = 1.5, y = 4.5, label = "True number of species = 5", size=5, color="blue")
#     p = p + theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold"))
# }
# plot.bpp.performance.pointrange = function(d, true.num.species) {
#     d$conversion.rate.level = factor(d$conversion.rate)
#     d$ymax = max(d$estimated.number.of.species)
#     p = ggplot(d, aes(x=conversion.rate.level, y=estimated.number.of.species))
#     p = p + stat_summary(fun.y = mean,
#                fun.ymin = function(a) min(a),
#                fun.ymax = function(a) max(a),
#                geom = "pointrange",
#                )
#     p = p + scale_y_continuous(limits=c(0, 20))
#     p = p + geom_hline(aes_string(yintercept=5), color="blue", linetype="dashed", size=0.75)
#     p = p + geom_text(data = NULL, x = 1.5, y = 4.5, label = "True number of species = 5", size=5, color="blue")
#     p = p + theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold"))
# }
# plot.bpp.performance.errorbar = function(d, true.num.species) {
#     d$conversion.rate.level = factor(d$conversion.rate)
#     d$ymax = max(d$estimated.number.of.species)
#     p = ggplot(d, aes(x=conversion.rate.level, y=estimated.number.of.species))
#     p = p + stat_summary(fun.y = mean,
#                fun.ymin = function(a) min(a),
#                fun.ymax = function(a) max(a),
#                geom = "errorbar",
#                width=0.2,
#                )
#     p = p + scale_y_continuous(limits=c(0, 20))
#     p = p + geom_hline(aes_string(yintercept=5), color="blue", linetype="dashed", size=0.75)
#     p = p + geom_text(data = NULL, x = 1.5, y = 4.5, label = "True number of species = 5", size=5, color="blue")
#     p = p + theme(axis.text=element_text(size=16), axis.title=element_text(size=14,face="bold"))
# }

# plot.all = function(d, true.num.species)  {
#     pdf("bpp-results-sp5.hist.pdf")
#     print(plot.bpp.performance.histograms(d, true.num.species))
#     print(dev.off())
#     print(pdf("bpp-results-sp5.boxplot.pdf"))
#     print(plot.bpp.performance.boxplot(d, true.num.species))
#     print(dev.off())
#     print(pdf("bpp-results-sp5.violin.pdf"))
#     print(plot.bpp.performance.violin(d, true.num.species))
#     print(dev.off())
#     print(pdf("bpp-results-sp5.scatter.pdf"))
#     print(plot.bpp.performance.scatter(d, true.num.species))
#     print(dev.off())
#     print(pdf("bpp-results-sp5.jitter.pdf"))
#     print(plot.bpp.performance.jitter(d, true.num.species))
#     print(dev.off())
#     print(pdf("bpp-results-sp5.pointrange.pdf"))
#     print(plot.bpp.performance.pointrange(d, true.num.species))
#     print(dev.off())
#     print(pdf("bpp-results-sp5.errorbar.pdf"))
#     print(plot.bpp.performance.errorbar(d, true.num.species))
#     dev.off()
# }

