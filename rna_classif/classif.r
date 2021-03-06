args = commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  args[1] = "C:\\Users\\Kristaps\\Documents\\GitHub\\RNAannotation\\rna_classif\\sorted.train.txt";
  args[2] = "C:\\Users\\Kristaps\\Documents\\GitHub\\RNAannotation\\rna_classif\\sorted.classes.txt";
  args[3] = "C:\\Users\\Kristaps\\Documents\\GitHub\\RNAannotation\\rna_classif\\classif_config.txt";
	#stop("Three arguments have to be passed (path to data, path to classes, path to config).", call.=FALSE)
}

if (!require("caret")) {
	  install.packages("caret")
}
library("caret")

if (!require("knitr")) {
	  install.packages("knitr")
}
library("knitr")

if (!require("markdown")) {
	  install.packages("markdown")
}
library("markdown")

# read the .csv files - it is not expected to have header data (data should start from the first line)
# the first column is expected to be the name of the entry - the rest are the values that will be used by the classifier
data <- read.table(args[1], sep="\t", header = FALSE ,row.names = 1)
# two columns expected - first one is the name of the entry, second one is the class of the entry
classes <- read.table(args[2], sep="\t", header = FALSE, row.names = 1)
# first column contains classfier name, second column is either empty or a parsable statement (data.frame with the tuning params for the classifier)
config <- read.table(args[3], sep="\t", header = FALSE)

# getting data for specific classes (need to figure out how to do this in a loop)
# ncData <- trainingSet[trainingSet$Class == 'ncRNA',]
# mData <- trainingSet[trainingSet$Class == 'mRNA',]
# otherData <- trainingSet[trainingSet$Class == 'other',]


# making the lists shorter (even out the class distribution)
# ncData <- ncData[1:300, ]
# mData <- mData[1:300, ]
# otherData <- otherData[1:300, ]

# hist(ncData$RFAM)
# hist(mData$RFAM)
# hist(otherData$RFAM)

# trainingSet <- rbind(ncData, mData, otherData)

# ncRNA <- ggplot(data[data$Expression < 2.5, ], aes(x = (Expression)))
# ncRNA + geom_density(aes(fill=factor(Class)), size=0.5)+facet_grid(Class ~ .) 

set.seed(100)

#if (length(args) < 4) {
	# it's doing something funky if the name column is left in (some sort of giant data structure gets built and I end up with 'Error: cannot allocate vector of size 46.0 Gb')
#	data$V1 <- NULL
#	data <- data.frame(data, Class = classes[, 2])
#	trainingClasses <- createDataPartition(data$Class, p = .75, list = FALSE)
#	trainingSet <- data[trainingClasses, ]
#	fullData <- data
#} else {
	# glueing the class column to the data
	data <- merge(x = data, y = classes, by = "row.names")
	data$Row.names <- NULL
	rownames(data) = data$Row.names
	colnames(data) = c('dir','ORFlength','SeqLength','GC','PFAM','Expression','RFAM','Class')

	# GC % instead of count
	data$GC <- with(data, GC / SeqLength)
	# ORF % instead of len
	data$ORFlength <- with(data, ORFlength / SeqLength)
	# PFAM and RFAM -1 -> 1 /// 10^-n -> log(of the val) (change the -1 to 1 before applying the log)
	data$PFAM[data$PFAM == -1] <- 1
	data$RFAM[data$RFAM == -1] <- 1
	data$PFAM <- with(data, log10(PFAM))
	data$RFAM <- with(data, log10(RFAM))
	
	fullData <- data
	
	# seqlength < 300
	data <- data[data$SeqLength < 300,]
	# a 1000 entries per attr (we're no longer interested in Dir)
	# order by orf % and pick them to be representative of the distribution - do the same for the rest of the cols
	gcData <- data
	gcData <- gcData[order(gcData$GC), ]
	gcData <- gcData[seq(1, nrow(gcData), floor(nrow(gcData) / 1000)), ]
	orfData <- data
	orfData <- orfData[order(orfData$ORFlength), ]
	orfData <- orfData[seq(1, nrow(orfData), floor(nrow(orfData) / 1000)), ]
	# Expression - representative of whole
	expData <- data
	expData <- expData[order(expData$Expression), ]
	expData <- expData[seq(1, nrow(expData), floor(nrow(expData) / 1000)), ]
	# PFAM and RFAM - drop 0 values - select a representative group from the stuff that's left
	pfamData <- data[data$PFAM != 0,]
	pfamData <- pfamData[order(pfamData$PFAM), ]
	pfamData <- pfamData[seq(1, nrow(pfamData), floor(nrow(pfamData) / 1000)), ]
	rfamData <- data[data$RFAM != 0,]
	rfamData <- rfamData[order(rfamData$RFAM), ]
	if (nrow(rfamData) > 1000) {
		rfamData <- rfamData[seq(1, nrow(rfamData), floor(nrow(rfamData) / 1000)), ]
	}

	trainingSet <- rbind(gcData, orfData, expData, pfamData, rfamData)
#}

fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10)

results <- list()

for (i in 1:nrow(config)) {
    set.seed(100)
    # checks if the second config column is empty and parses the grid function if it's there
    if (as.character(config[i, 2] != "")) {
        results[i] <- list(train(Class ~ ., data = trainingSet, 
                          method = as.character(config[i, 1]), 
                          trControl = fitControl,
                          verbose = FALSE,
                          tuneGrid = eval(parse(text=as.character(config[i, 2])))))

    } else {
        results[i] <- list(train(Class ~ ., data = trainingSet, 
                          method = as.character(config[i, 1]), 
                          trControl = fitControl,
                          verbose = FALSE))
    }
}

for (i in 1:nrow(config)) {
  cat("\n\nResults for ", config[i, 1], "\n\n")
  
  print(results[[i]]) 
  x11()
  plot(results[[i]]) 
  
  summary(results[[i]])
  
  # classifying everything we have
  x_test <- data[, 1:7]
  y_test <- data[, 8]
  predictions <- predict(results[[i]], x_test)
  print(confusionMatrix(predictions, y_test))
  
  print(results[[i]]$finalModel)
}
#knitr::opts_chunk$set(echo=FALSE, fig.path='chunks/plot-', cache=TRUE)
#knit("classif_report.rmd")
#markdownToHTML('classif_report.md', 'classif_report.html', options=c("use_xhml"))