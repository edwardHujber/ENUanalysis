library("ggplot2")

# filtereddat1_filtered <- filtereddat1[which(filtereddat1$altcount>8),]
# dat_new_filtered <- dat_new[which(dat_new$altcount>8 & dat_new$HR<3.5),]

cols <- c("filtereddat1"="green","filtereddat2"="orange","filtered"="green", "new_filtered"="orange")

# flag < 2.5  more stringent than recomended
ggplot(filtereddat1,aes(as.numeric(QD),colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   geom_vline(xintercept= 2 ) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   ggtitle("QD")
   

# flag > 30  more stringent than recomended
ggplot(filtereddat1,aes(FS,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   coord_cartesian( ylim = c(0, 0.001))+
   geom_vline(xintercept= 30 ) +
   ggtitle("FS, y-zoomed")

ggplot(filtereddat1,aes(FS,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   coord_cartesian( xlim = c(0, 30))+
   geom_vline(xintercept= 30 ) +
   ggtitle("FS, x-zoomed")


# flag > 3  as stringent as recomended
ggplot(filtereddat1,aes(SOR,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   geom_vline(xintercept= 3 ) +
   ggtitle("SOR")


# flag < 55  more stringent than recomended
ggplot(filtereddat1,aes(MQ,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   coord_cartesian( xlim = c(50, 60))+
   geom_vline(xintercept= 40 ) +
   ggtitle("MQ, x-zoomed")

ggplot(filtereddat1,aes(MQ,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   coord_cartesian( ylim = c(0, 0.01))+
   geom_vline(xintercept= 40 ) +
   ggtitle("MQ, y-zoomed")

# flag < -2.5 or > 2.5  as stringent as recomended
ggplot(filtereddat1,aes(as.numeric(MQRankSum),colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   geom_vline(xintercept= -2.5 ) +
   geom_vline(xintercept= 2.5 ) +
   ggtitle("MQRankSum")


# flag < -2 or > 2  more stringent than recomended
ggplot(filtereddat1,aes(as.numeric(ReadPosRankSum),colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   geom_vline(xintercept= -2 ) +
   geom_vline(xintercept= 2 ) +
   ggtitle("ReadPosRankSum")

################################################################


ggplot(filtereddat1,aes(VAF,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
  # geom_density(data=filtereddat, aes(colour="filtered", fill="filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   geom_vline(xintercept= 3.5 ) +
   coord_cartesian( xlim = c(0, 1))+
   ggtitle("VAF")

################################################################

ggplot(filtereddat1,aes(DP,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   coord_cartesian( xlim = c(0, 30))+
   ggtitle("DP")

ggplot(filtereddat1,aes(altcount,colour="filtereddat1", fill="filtereddat1")) + 
   geom_density(bw = 1.0, alpha=0.2, size=1.5) + 
   geom_density(data=filtereddat2, bw=1.0, aes(colour="filtereddat2", fill="filtereddat2"),alpha=0.2, size=1.5) +
   # geom_density(data=filtereddat1_filtered, aes(colour="old_filtered", fill="old_filtered"),alpha=0.2, size=1.5) +
   # geom_density(data=dat_new_filtered, aes(colour="new_filtered", fill="new_filtered"),alpha=0.2, size=1.5) +
   scale_colour_manual(name="Datasets",values=cols) +
   scale_fill_manual(name="Datasets",values=cols) +
   coord_cartesian( xlim = c(0, 15))+
   ggtitle("altcount")


ggplot(filtereddat1, aes(HR, altcount)) +
   geom_point(alpha=0.5) +
   geom_point(data=filtereddat2, color="red",alpha=0.5) +
   scale_x_log10()  +
   scale_y_log10() +   
   geom_vline(xintercept= 3.5 ) +   
   geom_hline(yintercept= 8)



