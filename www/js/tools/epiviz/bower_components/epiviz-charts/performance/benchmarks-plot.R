scaled <- jsonlite::fromJSON("~/Desktop/projects/epiviz/epiviz-chart/performance/results/scaled.log")

datapoints <- names(scaled)

shttpTime <- sapply(names(scaled), function(dp) {
  mean(scaled[[dp]][['total_http_time']])
})

shttpTime_sd <- sapply(names(scaled), function(dp) {
  sd(scaled[[dp]][['total_http_time']])
})

sdrawTime <- sapply(names(scaled), function(dp) {
  mean(scaled[[dp]][['total_draw_time']])
})

sdrawTime_sd <- sapply(names(scaled), function(dp) {
  sd(scaled[[dp]][['total_draw_time']])
})

slatencyTime <- sapply(names(scaled), function(dp) {
  mean(scaled[[dp]][['total_latency_time']])
})

slatencyTime_sd <- sapply(names(scaled), function(dp) {
  sd(scaled[[dp]][['total_latency_time']])
})

srequestSize <- sapply(names(scaled), function(dp) {
  mean(scaled[[dp]][['total_request_size']])/1000
})

srequestSize_sd <- sapply(names(scaled), function(dp) {
  sd(scaled[[dp]][['total_request_size']])
})

unscaled <- jsonlite::fromJSON("~/Desktop/projects/epiviz/epiviz-chart/performance/results/unscaled.log")
datapoints <- names(unscaled)

uhttpTime <- sapply(names(unscaled), function(dp) {
  mean(unscaled[[dp]][['total_http_time']])
})

uhttpTime_sd <- sapply(names(unscaled), function(dp) {
  sd(unscaled[[dp]][['total_http_time']])
})

udrawTime <- sapply(names(unscaled), function(dp) {
  mean(unscaled[[dp]][['total_draw_time']])
})

udrawTime_sd <- sapply(names(unscaled), function(dp) {
  sd(unscaled[[dp]][['total_draw_time']])
})

ulatencyTime <- sapply(names(unscaled), function(dp) {
  mean(unscaled[[dp]][['total_latency_time']])
})

ulatencyTime_sd <- sapply(names(unscaled), function(dp) {
  sd(unscaled[[dp]][['total_latency_time']])
})

urequestSize <- sapply(names(unscaled), function(dp) {
  mean(unscaled[[dp]][['total_request_size']]) /1000
})

urequestSize_sd <- sapply(names(unscaled), function(dp) {
  sd(unscaled[[dp]][['total_request_size']])
})

formatted <- data.frame(grange = datapoints,
                        shttpTime = shttpTime,
                        shttpTime_sd = shttpTime_sd,
                        sdrawTime = sdrawTime,
                        sdrawTime_sd = sdrawTime_sd,
                        slatencyTime = slatencyTime,
                        slatencyTime_sd = slatencyTime_sd,
                        srequestSize = srequestSize,
                        uhttpTime = uhttpTime,
                        uhttpTime_sd = uhttpTime_sd,
                        udrawTime = udrawTime,
                        udrawTime_sd = udrawTime_sd,
                        ulatencyTime = ulatencyTime,
                        ulatencyTime_sd = ulatencyTime_sd,
                        urequestSize = urequestSize,
                        stringsAsFactors = FALSE)

library(ggplot2)
library(gridExtra)

drawTimePlot <- ggplot(data=formatted, aes(x = grange, group = 2)) +
  geom_line(aes(y = udrawTime, colour = "unsummarized")) +
  geom_line(aes(y = sdrawTime, colour = "summarized")) +
  ylab("draw time (ms)") +
  xlab("genomic range (bp)") +
  geom_errorbar(aes(ymin=udrawTime-udrawTime_sd, ymax=udrawTime+udrawTime_sd, colour = "unsummarized"), width=.2) +
  geom_errorbar(aes(ymin=sdrawTime-sdrawTime_sd, ymax=sdrawTime+sdrawTime_sd, colour = "summarized"), width=.2) +
  theme(legend.position="top") + 
  ggtitle(" ") +
  scale_x_discrete(limits= formatted$grange)

library(reshape2)
library(data.table)
barDatahttp <- melt(formatted[, c("grange", "shttpTime", "uhttpTime")], id.vars=1)
barDatalatency <- melt(formatted[, c("grange", "slatencyTime", "ulatencyTime")], id.vars=1)
ftable <- as.data.table(formatted)

ftable$sdataTime = ftable$shttpTime - ftable$slatencyTime
ftable$udataTime = ftable$uhttpTime - ftable$ulatencyTime

barData <- melt(ftable[, c("grange", "urequestSize", "srequestSize", "shttpTime_sd",
                           "uhttpTime_sd", "slatencyTime", "ulatencyTime", "sdataTime", 
                           "udataTime", "shttpTime", "uhttpTime")],
                measure = list(c("shttpTime", "uhttpTime"),
                               c("sdataTime", "udataTime"),
                               c("slatencyTime", "ulatencyTime"),
                               c("shttpTime_sd", "uhttpTime_sd"),
                               c("srequestSize", "urequestSize")),
                value.name = c("httpTime", "dataTime", "latencyTime", "httpTime_sd", "requestSize"))

barData[barData$variable == 1]$variable <- "summarized"
barData[barData$variable == 2 ]$variable <- "unsummarized"
bData <- melt(barData, id=c("grange", "variable"))
names(bData) <- c("grange", "type", "benchmark", "time")
bData[bData$grange == "10K" ]$grange <- 100
bData[bData$grange == "100K" ]$grange <- 150
bData[bData$grange == "1M" ]$grange <- 200
bData[bData$grange == "10M" ]$grange <- 250
bData[bData$grange == "100M" ]$grange <- 300
bData[bData$grange == "chr" ]$grange <- 350
bData$grange <- as.numeric(bData$grange)

sumdataTime <- bData[bData$type == "summarized" & bData$benchmark %in% c("dataTime", "latencyTime")]
usumdataTime <- bData[bData$type == "unsummarized" & bData$benchmark %in% c("dataTime", "latencyTime")]
sumhttpTime <- bData[bData$type == "summarized" & bData$benchmark == "httpTime"]
usumhttpTime <- bData[bData$type == "unsummarized" & bData$benchmark == "httpTime"]

stackedBPlot <-  ggplot(group = 2) +
  geom_col(data = sumdataTime,
           aes(x = grange - 10.5, y = time/1000, fill=benchmark),
           width = 20, colour="grey40", alpha = 0.9) + 
  geom_col(data = usumdataTime,
           aes(x = grange + 10.5, y = time/1000, fill=benchmark),
           width = 20, colour="grey40", alpha = 0.9) + 
  geom_text(data=bData[bData$type == "summarized" & bData$benchmark == "requestSize"],
            aes(label=paste0(time, " KB"), x = grange - 10.5, y = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5)),
            size=3.5, angle = 90) +
  geom_text(data=bData[bData$type == "unsummarized" & bData$benchmark == "requestSize"],
            aes(label=paste0(time, " KB"), x = grange + 10.5, y = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5)),
            size=3.5, angle = 90) +
  geom_text(data=bData[bData$type == "summarized" & bData$benchmark == "requestSize"],
            aes(label="s", x = grange - 10.5, y = c(-0.2, -0.2, -0.2, -0.2, -0.2, -0.2)),
            size=3.5) +
  geom_text(data=bData[bData$type == "unsummarized" & bData$benchmark == "requestSize"],
            aes(label="u", x = grange + 10.5, y = c(-0.2, -0.2, -0.2, -0.2, -0.2, -0.2)),
            size=3.5) +
  geom_errorbar(data=bData[bData$type == "summarized" & bData$benchmark == "httpTime_sd"],
                aes(x = grange - 10.5, ymin = (sumhttpTime$time - time)/1000, 
                    ymax = (sumhttpTime$time + time)/1000), width=10,
                position=position_dodge(.9)) +
  geom_errorbar(data=bData[bData$type == "unsummarized" & bData$benchmark == "httpTime_sd"],
                aes(x = grange + 10.5, ymin = (usumhttpTime$time - time)/1000, 
                    ymax = (usumhttpTime$time + time)/1000), width=10,
                position=position_dodge(.9)) +
  theme(legend.position="top", plot.title = element_text(hjust = 0.5, size=10)) + 
  ggtitle("s = summarized, u = unsummarized") +
  ylab("mean http time (s)") +
  scale_x_continuous(breaks=c(100, 150, 200, 250, 300, 350), name = "genomic range (bp)", labels = c("10K", "100K",
                                               "1M", "10M", "100M", "chr"))

grid.arrange(drawTimePlot, stackedBPlot, nrow = 1)

