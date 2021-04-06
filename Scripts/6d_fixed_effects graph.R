require(ggplot2)

alpha2
alpha2 <- alpha2[c(1,7,6,10,9,8,5,4,3,2,11),]
alpha2$type <- factor(alpha2$type, levels =  c("primary","secondary","pasture","crop","urban"))
alpha2$intensity <- factor(alpha2$intensity, levels = c("Minimal","Light", "Intense", "All"))

pd <- position_dodge(0.5)## this is to make the point on one land use separated
ggplot(alpha2, aes(x=type, y=Gow.Rao, colour=intensity, group=intensity)) + 
  geom_errorbar(aes(ymin=Gow.Rao + (lower.CI - Value) , ymax=Gow.Rao+ (upper.CI - Value)), colour="black", width=.1, position=pd, linetype = 1) +
  geom_point(position=pd,size=6)+
  xlab("Land use") +
  ylab("Estimate Rao'Q") +
  scale_colour_hue(name="Use intensity",    # Legend label, use darker colors
                   breaks=c("Minimal", "Light","Intense","All"),
                   labels=c("Minimal", "Light","Intense","All"),
                   l=40) +                    # Use darker colors, lightness=40+
  expand_limits(y=0) +
  ylim(0.275,0.40) +
  theme_classic() +
  theme(legend.justification=c(1,0),
        legend.position=c(0.95,0.65))


alpha2$Gow.Rao + (alpha2$lower.CI - alpha2$Value)
