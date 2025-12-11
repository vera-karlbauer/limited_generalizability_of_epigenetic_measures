### Title: "Cross-tissue correlations: Define utilities"
### Author: "Vera N. Karlbauer"
### Contact: "vera_karlbauer@psych.mpg.de"
### Date created: "2025-11-25"
### Purpose: Define utlities for cross-tissue project

### Define color palettes based on color brewer

## grey-to-purple palette (for cross-tissue data)
# get RdGy and PRGn palettes
rdgy_colors <- brewer.pal(11, "RdGy")
prgn_colors <- brewer.pal(11, "PRGn")
# combine right side of RdGy and left side of PRGn
GyPr_palette <- colorRampPalette(c(rev(rdgy_colors[6:11]), rev(prgn_colors[1:5])))(200)

## grey-to-red palette (for blood-only data)
GyRd_palette <- colorRampPalette(rev(rdgy_colors))(200)

## grey-to-blue palette (for saliva-only data)
# get RdBu palette
rdbu_colors <- brewer.pal(11, "RdBu")
# combine right side of RdGy and right side of RdBu
GyBu_palette <- colorRampPalette(c(rev(rdgy_colors[6:11]), rdbu_colors[7:11]))(200)
# remove intermediary palettes
remove(rdgy_colors, prgn_colors, rdbu_colors)


### Define order of epigenetic clocks for plotting
clockorder = c("Horvath", "PCHorvath", "SkinBlood", "PCSkinBlood", "Hannum", "PCHannum", "PedBE", "Wu", "PhenoAge", "PCPhenoAge", "GrimAge", "PCGrimAge", "GrimAge2", "DunedinPACE")

### Define order of epigenetic scores for plotting


