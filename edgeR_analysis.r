setwd("/Users/edwards/RICE/Lund+Soil/greenhouse/edgeR/")
load("comp.site.whole.rda")
library(gplots)
tax <- read.table(file = "../../q.tax", header = T, row.names = 1)
tax$OTU <- paste("Otu", row.names(tax), sep = "")
row.names(tax) <- NULL


da.counts <- data.frame(Value = c(
nrow(subset(comp.site.whole, Comp == "Endosphere" & Site == "Arbuckle" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizoplane" & Site == "Arbuckle" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizosphere" & Site == "Arbuckle" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Endosphere" & Site == "Arbuckle" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizoplane" & Site == "Arbuckle" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizosphere" & Site == "Arbuckle" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Endosphere" & Site == "Davis" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizoplane" & Site == "Davis" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizosphere" & Site == "Davis" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Endosphere" & Site == "Davis" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizoplane" & Site == "Davis" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizosphere" & Site == "Davis" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Endosphere" & Site == "Sacramento" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizoplane" & Site == "Sacramento" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizosphere" & Site == "Sacramento" & logFC > 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Endosphere" & Site == "Sacramento" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizoplane" & Site == "Sacramento" & logFC < 0 & padj <= 0.05)),
nrow(subset(comp.site.whole, Comp == "Rhizosphere" & Site == "Sacramento" & logFC < 0 & padj <= 0.05))),
Site = c(rep("Arbuckle", 6), rep("Davis", 6), rep("Sacramento", 6)),
Comp= rep(c("Endosphere", "Rhizoplane", "Rhizosphere"), 6),
Category = rep(c(rep("up" , 3), rep("down", 3)), 3))

da.counts$x.vp <- rep(c(rep(7.5 , 3), rep(-10, 3)), 3)
da.counts$y.vp <- rep(110, nrow(da.counts))
da.counts$x.ma <- rep(1, nrow(da.counts))
da.counts$y.ma <- rep(c(rep(10 , 3), rep(-10, 3)), 3)

ggplot(comp.site.whole, aes(x = logCPM, y = logFC, color = ifelse(padj > 0.05, "ns", ifelse(logFC < 0, "up", "down")))) +
	geom_point(alpha = 1/3, size = 1.5) +
	scale_color_manual(values = c("red", "grey32", "blue"), guide = FALSE) +
	facet_grid(Site ~ Comp) +
	geom_text(aes(x = x.ma, y = y.ma, label = Value, group= NULL), data= da.counts, inherit.aes = FALSE, parse = FALSE, size= 10, color = c('red', 'blue')) +
	theme_bw() +
	theme(text = element_text(size = 20))

ggplot(comp.site.whole, aes(x = logFC, y = -log10(padj), color = ifelse(padj > 0.05, "ns", ifelse(logFC < 0, "up", "down")))) +
	geom_vline(x = 0, color = 'black') +
	geom_point(alpha = 1/3, size = 1.5) +
	theme_bw() +
	scale_color_manual(values = c("red", "grey32", "blue"), guide = FALSE) +
	facet_grid(Site ~ Comp) +
	geom_text(aes(x = x.vp, y = y.vp, label = Value, group= NULL), data= da.counts, inherit.aes = FALSE, parse = FALSE, size = 8, color= c('red', 'blue')) +
	labs(x = "Log10 Fold Change", y = "-Log10 Adjusted P Value") +
	theme(text = element_text(size = 20))

comp.site.whole$Mean2 <- (2 * (10^comp.site.whole$logCPM)) / (1 + (10^comp.site.whole$logFC))
comp.site.whole$Mean1 <- comp.site.whole$Mean2 * (10^comp.site.whole$logFC)

ggplot(comp.site.whole, aes(x = log10(Mean1), y = log10(Mean2), color = ifelse(padj >= 0.05, "ns", "s"))) +
	geom_point(alpha = 1/3) +
	theme_bw() +
	scale_color_manual(values = c("grey32", "red"), guide = FALSE) +
	facet_grid(Site ~ Comp)

arb <- subset(comp.site.whole, padj <= 0.05 & Site == "Arbuckle" & logFC > 0)
arb <- arb[!duplicated(arb$OTU),]
arb.e.up <- subset(arb, Comp == "Endosphere")
arb.rp.up <- subset(arb, Comp == "Rhizoplane")
arb.rs.up <- subset(arb, Comp == "Rhizosphere")
arb.down <- subset(comp.site.whole, padj <= 0.05 & Site == "Arbuckle" & logFC < 0)
arb.down <- arb.down[!duplicated(arb.down$OTU),]
arb.e.down <- subset(arb.down, Comp == "Endosphere")
arb.rp.down <- subset(arb.down, Comp == "Rhizoplane")
arb.rs.down <- subset(arb.down, Comp == "Rhizosphere")

dav <- subset(comp.site.whole, padj <= 0.05 & Site == "Davis" & logFC > 0)
dav <- dav[!duplicated(dav$OTU),]
dav.e.up <- subset(dav, Comp == "Endosphere")
dav.rp.up <- subset(dav, Comp == "Rhizoplane")
dav.rs.up <- subset(dav, Comp == "Rhizosphere")
dav.down <- subset(comp.site.whole, padj <= 0.05 & Site == "Davis" & logFC < 0)
dav.down <- dav.down[!duplicated(dav.down$OTU),]
dav.e.down <- subset(dav.down, Comp == "Endosphere")
dav.rp.down <- subset(dav.down, Comp == "Rhizoplane")
dav.rs.down <- subset(dav.down, Comp == "Rhizosphere")

sac <- subset(comp.site.whole, padj <= 0.05 & Site == "Sacramento" & logFC > 0)
sac <- sac[!duplicated(sac$OTU),]
sac.e.up <- subset(sac, Comp == "Endosphere")
sac.rp.up <- subset(sac, Comp == "Rhizoplane")
sac.rs.up <- subset(sac, Comp == "Rhizosphere")
sac.down <- subset(comp.site.whole, padj <= 0.05 & Site == "Sacramento" & logFC < 0)
sac.down <- sac.down[!duplicated(sac.down$OTU),]
sac.e.down <- subset(sac.down, Comp == "Endosphere")
sac.rp.down <- subset(sac.down, Comp == "Rhizoplane")
sac.rs.down <- subset(sac.down, Comp == "Rhizosphere")

e.up <- list(Arbuckle = arb.e.up$OTU, Davis = dav.e.up$OTU, Sacramento= sac.e.up$OTU)
e.down <-  list(Arbuckle = arb.e.down$OTU, Davis = dav.e.down$OTU, Sacramento= sac.e.down$OTU)
venn.diagram(e.up, fill = c("brown4", "darkorange", "grey50"), filename = "e.up.venn.tiff", cex= 3)
venn.diagram(e.down, fill = c("brown4", "darkorange", "grey50"), filename = "e.down.venn.tiff", cex= 3)

rp.up <- list(Arbuckle = arb.rp.up$OTU, Davis = dav.rp.up$OTU, Sacramento= sac.rp.up$OTU)
rp.down <-  list(Arbuckle = arb.rp.down$OTU, Davis = dav.rp.down$OTU, Sacramento= sac.rp.down$OTU)
venn.diagram(rp.up, fill = c("brown4", "darkorange", "grey50"), filename = "rp.up.venn.tiff", cex= 3)
venn.diagram(rp.down, fill = c("brown4", "darkorange", "grey50"), filename = "rp.down.venn.tiff", cex= 2)

rs.up <- list(Arbuckle = arb.rs.up$OTU, Davis = dav.rs.up$OTU, Sacramento= sac.rs.up$OTU)
rs.down <- list(Arbuckle = arb.rs.down$OTU, Davis = dav.rs.down$OTU, Sacramento= sac.rs.down$OTU)
venn.diagram(rs.up, fill = c("brown4", "darkorange", "grey50"), filename = "rs.up.venn.tiff", cex= 3)
venn.diagram(rs.down, fill = c("brown4", "darkorange", "grey50"), filename = "rs.down.venn.tiff",cex= 3)

core.e <- intersect(arb.e.up$OTU, intersect(dav.e.up$OTU, sac.e.up$OTU))
core.e.tax <- tax[match(core.e, tax$OTU),]
ggplot(core.e.tax, aes(x= factor(""), fill = Phylum)) +
	geom_bar(width = 1) +
	coord_polar(theta = "y")
	
arb.e.tax <- cbind(tax[match(arb.e.up$OTU, tax$OTU),], Site = rep("Arbuckle", nrow(arb.e.up)), Compartment = rep("Endosphere", nrow(arb.e.up)))
dav.e.tax <- cbind(tax[match(dav.e.up$OTU, tax$OTU),], Site = rep("Davis", nrow(dav.e.up)), Compartment = rep("Endosphere", nrow(dav.e.up)))
sac.e.tax <- cbind(tax[match(sac.e.up$OTU, tax$OTU),], Site = rep("Sacramento", nrow(sac.e.up)), Compartment = rep("Endosphere", nrow(sac.e.up)))

arb.rp.tax <- cbind(tax[match(arb.rp.up$OTU, tax$OTU),], Site = rep("Arbuckle", nrow(arb.rp.up)), Compartment = rep("Rhizoplane", nrow(arb.rp.up)))
dav.rp.tax <- cbind(tax[match(dav.rp.up$OTU, tax$OTU),], Site = rep("Davis", nrow(dav.rp.up)), Compartment = rep("Rhizoplane", nrow(dav.rp.up)))
sac.rp.tax <- cbind(tax[match(sac.rp.up$OTU, tax$OTU),], Site = rep("Sacramento", nrow(sac.rp.up)), Compartment = rep("Rhizoplane", nrow(sac.rp.up)))

arb.rs.tax <- cbind(tax[match(arb.rs.up$OTU, tax$OTU),], Site = rep("Arbuckle", nrow(arb.rs.up)), Compartment = rep("Rhizosphere", nrow(arb.rs.up)))
dav.rs.tax <- cbind(tax[match(dav.rs.up$OTU, tax$OTU),], Site = rep("Davis", nrow(dav.rs.up)), Compartment = rep("Rhizosphere", nrow(dav.rs.up)))
sac.rs.tax <- cbind(tax[match(sac.rs.up$OTU, tax$OTU),], Site = rep("Sacramento", nrow(sac.rs.up)), Compartment = rep("Rhizosphere", nrow(sac.rs.up)))

whole.tax <- rbind(arb.e.tax, dav.e.tax, sac.e.tax, arb.rp.tax, dav.rp.tax, sac.rp.tax, arb.rs.tax, dav.rs.tax, sac.rs.tax)
top_15 <- names(sort(table(whole.tax$Phylum), decreasing = T))[1:15]

arb.e.tax2 <- cbind(as.data.frame(table(arb.e.tax$Phylum)), Compartment = "Endosphere", Site = "Arbuckle")
arb.rp.tax2 <- cbind(as.data.frame(table(arb.rp.tax$Phylum)), Compartment = "Rhizoplane", Site = "Arbuckle")
arb.rs.tax2 <- cbind(as.data.frame(table(arb.rs.tax$Phylum)), Compartment = "Rhizosphere", Site = "Arbuckle")
dav.e.tax2 <- cbind(as.data.frame(table(dav.e.tax$Phylum)), Compartment = "Endosphere", Site = "Davis")
dav.rp.tax2 <- cbind(as.data.frame(table(dav.rp.tax$Phylum)), Compartment = "Rhizoplane", Site = "Davis")
dav.rs.tax2 <- cbind(as.data.frame(table(dav.rs.tax$Phylum)), Compartment = "Rhizosphere", Site = "Davis")
sac.e.tax2 <- cbind(as.data.frame(table(sac.e.tax$Phylum)), Compartment = "Endosphere", Site = "Sacramento")
sac.rp.tax2 <- cbind(as.data.frame(table(sac.rp.tax$Phylum)), Compartment = "Rhizoplane", Site = "Sacramento")
sac.rs.tax2 <- cbind(as.data.frame(table(sac.rs.tax$Phylum)), Compartment = "Rhizosphere", Site = "Sacramento")

whole.tax2 <- rbind(arb.e.tax2, dav.e.tax2, sac.e.tax2, arb.rp.tax2, dav.rp.tax2, sac.rp.tax2, arb.rs.tax2, dav.rs.tax2, sac.rs.tax2)
names(whole.tax2)[1] <- "Phylum"
top.15.names <- names(sort(tapply(whole.tax2$Freq, whole.tax2$Phylum, sum), decreasing = T))[1:15]
top.15 <- whole.tax2[whole.tax2$Phylum%in%top.15.names,]
others <- whole.tax2[!whole.tax2$Phylum%in%top.15.names,]
other <- data.frame(Phylum = rep("Other", 9), Freq = as.vector(tapply(others$Freq, list(others$Compartment, others$Site), sum)), Compartment = rep(c("Endosphere", "Rhizoplane", "Rhizosphere"), 3), Site = c(rep("Arbuckle", 3), rep("Davis", 3), rep("Sacramento", 3)))
all.tax2 <- rbind(top.15, other)

plot.cols <- c(brewer.pal(8, "Set1"), brewer.pal(7, "Set2"), "grey50")
ggplot(all.tax2, aes(x = factor(""), y = Freq, fill = Phylum)) + 
	geom_bar(width = 1, position = "fill", stat = "identity") +
	coord_polar(theta = "y") +
	scale_fill_manual(values = plot.cols) +
	facet_grid(Site~Compartment) +
	theme_bw() +
	theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
	labs(x = "", y = "")
	
arb.e.tax.down <- cbind(tax[match(arb.e.down$OTU, tax$OTU),], Site = rep("Arbuckle", nrow(arb.e.down)), Compartment = rep("Endosphere", nrow(arb.e.down)))
dav.e.tax.down <- cbind(tax[match(dav.e.down$OTU, tax$OTU),], Site = rep("Davis", nrow(dav.e.down)), Compartment = rep("Endosphere", nrow(dav.e.down)))
sac.e.tax.down <- cbind(tax[match(sac.e.down$OTU, tax$OTU),], Site = rep("Sacramento", nrow(sac.e.down)), Compartment = rep("Endosphere", nrow(sac.e.down)))

arb.rp.tax.down <- cbind(tax[match(arb.rp.down$OTU, tax$OTU),], Site = rep("Arbuckle", nrow(arb.rp.down)), Compartment = rep("Rhizoplane", nrow(arb.rp.down)))
dav.rp.tax.down <- cbind(tax[match(dav.rp.down$OTU, tax$OTU),], Site = rep("Davis", nrow(dav.rp.down)), Compartment = rep("Rhizoplane", nrow(dav.rp.down)))
sac.rp.tax.down <- cbind(tax[match(sac.rp.down$OTU, tax$OTU),], Site = rep("Sacramento", nrow(sac.rp.down)), Compartment = rep("Rhizoplane", nrow(sac.rp.down)))

arb.rs.tax.down <- cbind(tax[match(arb.rs.down$OTU, tax$OTU),], Site = rep("Arbuckle", nrow(arb.rs.down)), Compartment = rep("Rhizosphere", nrow(arb.rs.down)))
dav.rs.tax.down <- cbind(tax[match(dav.rs.down$OTU, tax$OTU),], Site = rep("Davis", nrow(dav.rs.down)), Compartment = rep("Rhizosphere", nrow(dav.rs.down)))
sac.rs.tax.down <- cbind(tax[match(sac.rs.down$OTU, tax$OTU),], Site = rep("Sacramento", nrow(sac.rs.down)), Compartment = rep("Rhizosphere", nrow(sac.rs.down)))

whole.tax.down <- rbind(arb.e.tax.down, dav.e.tax.down, sac.e.tax.down, arb.rp.tax.down, dav.rp.tax.down, sac.rp.tax.down, arb.rs.tax.down, dav.rs.tax.down, sac.rs.tax.down)
top_15.down <- names(sort(table(whole.tax.down$Phylum), decreasing = T))[1:15]

arb.e.tax2.down <- cbind(as.data.frame(table(arb.e.tax.down$Phylum)), Compartment = "Endosphere", Site = "Arbuckle")
arb.rp.tax2.down <- cbind(as.data.frame(table(arb.rp.tax.down$Phylum)), Compartment = "Rhizoplane", Site = "Arbuckle")
arb.rs.tax2.down <- cbind(as.data.frame(table(arb.rs.tax.down$Phylum)), Compartment = "Rhizosphere", Site = "Arbuckle")
dav.e.tax2.down <- cbind(as.data.frame(table(dav.e.tax.down$Phylum)), Compartment = "Endosphere", Site = "Davis")
dav.rp.tax2.down <- cbind(as.data.frame(table(dav.rp.tax.down$Phylum)), Compartment = "Rhizoplane", Site = "Davis")
dav.rs.tax2.down <- cbind(as.data.frame(table(dav.rs.tax.down$Phylum)), Compartment = "Rhizosphere", Site = "Davis")
sac.e.tax2.down <- cbind(as.data.frame(table(sac.e.tax.down$Phylum)), Compartment = "Endosphere", Site = "Sacramento")
sac.rp.tax2.down <- cbind(as.data.frame(table(sac.rp.tax.down$Phylum)), Compartment = "Rhizoplane", Site = "Sacramento")
sac.rs.tax2.down <- cbind(as.data.frame(table(sac.rs.tax.down$Phylum)), Compartment = "Rhizosphere", Site = "Sacramento")

whole.tax2.down <- rbind(arb.e.tax2.down, dav.e.tax2.down, sac.e.tax2.down, arb.rp.tax2.down, dav.rp.tax2.down, sac.rp.tax2.down, arb.rs.tax2.down, dav.rs.tax2.down, sac.rs.tax2.down)
names(whole.tax2.down)[1] <- "Phylum"
top.15.names.down <- names(sort(tapply(whole.tax2.down$Freq, whole.tax2.down$Phylum, sum), decreasing = T))[1:15]
top.15.down <- whole.tax2.down[whole.tax2.down$Phylum%in%top.15.names.down,]
others.down <- whole.tax2.down[!whole.tax2.down$Phylum%in%top.15.names.down,]
other.down <- data.frame(Phylum = rep("Other", 9), Freq = as.vector(tapply(others.down$Freq, list(others.down$Compartment, others.down$Site), sum)), Compartment = rep(c("Endosphere", "Rhizoplane", "Rhizosphere"), 3), Site = c(rep("Arbuckle", 3), rep("Davis", 3), rep("Sacramento", 3)))
all.tax2.down <- rbind(top.15.down, other.down)

plot.cols <- c(brewer.pal(8, "Set1"), brewer.pal(7, "Set2"), "grey50")
ggplot(all.tax2.down, aes(x = factor(""), y = Freq, fill = Phylum)) + 
	geom_bar(width = 1, position = "fill", stat = "identity") +
	coord_polar(theta = "y") +
	scale_fill_manual(values = plot.cols) +
	facet_grid(Site~Compartment) +
	theme_bw() +
	theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20)) +
	labs(x = "", y = "")


