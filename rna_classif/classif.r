args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  # adding hard-coded control args to make it easier to run form the ide
  args[1] = "C:\\Users\\Kristaps\\Desktop\\sorted.train.txt";
  args[2] = "C:\\Users\\Kristaps\\Desktop\\sorted.classes.txt";
  args[3] = "C:\\Users\\Kristaps\\Desktop\\classif_config.txt";
	#stop("Three arguments have to be passed (path to data, path to classes, path to config).", call.=FALSE)
}


if (!require("caret")) {
	install.packages("caret")
}
library("caret")

# read the .csv files - it is not expected to have header data (data should start from the first line)
# the first column is expected to be the name of the entry - the rest are the values that will be used by the classifier
data <- read.table(args[1], sep="\t", header = FALSE)
# two columns expected - first one is the name of the entry, second one is the class of the entry
classes <- read.table(args[2], sep="\t", header = FALSE)
# first column contains classfier name, second column is either empty or a parsable statement (data.frame with the tuning params for the classifier)
config <- read.table(args[3], sep="\t", header = FALSE)

# it's doing something funky if the name column is left in (some sort of giant data structure gets built and I end up with 'Error: cannot allocate vector of size 46.0 Gb')
data$V1 <- NULL

# glueing the class column to the data
data <- data.frame(data, Class = classes[, 2])

fullData <- data

# ORF < 200
data <- data[data$V3 < 200,]

# getting data for specific classes (need to figure out how to do this in a loop)
ncData <- data[data$Class == 'ncRNA',]
mData <- data[data$Class == 'mRNA',]
otherData <- data[data$Class == 'other',]

# making the lists shorter (even out the class distribution)
ncData <- ncData[1:100, ]
mData <- mData[1:100, ]
otherData <- otherData[1:100, ]

data <- rbind(ncData, mData, otherData)

set.seed(100)

trainingClasses <- createDataPartition(data$Class, p = .75, list = FALSE)
trainingSet <- data[trainingClasses, ]
testSet <- data[-trainingClasses, ]
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
    print(results[[i]]) 
    plot(results[[i]]) 
    
    # info about the classifier
    # gbm training with 225 entries (~50% accuracy)
    # V7  V7 39.957349 (Exp_level)
    # V3  V3 27.171016 (ORFlength)
    # V4  V4 18.636434 (SeqLength)
    # V5  V5 13.225915 (GC_count)
    # V2  V2  1.009286 (Dir)
    # V6  V6  0.000000 (PFAM_E-value)
    # V8  V8  0.000000 (RFAM_E-Value)
    summary(results[[i]])
    
    # classifying the test set
    x_test <- testSet[, 1:7]
    y_test <- testSet[, 8]
    predictions <- predict(results[[i]], x_test)
    confusionMatrix(predictions, y_test)
    
    # classifying everything we have
    x_test <- fullData[, 1:7]
    y_test <- fullData[, 8]
    predictions <- predict(results[[i]], x_test)
    confusionMatrix(predictions, y_test)
    
    # Reference
    # Prediction  mRNA ncRNA other
    # mRNA  23031  3588 21083
    # ncRNA  2782  1311  7791
    # other   996  1195 16779
    #
    # Overall Statistics
    #
    # Accuracy : 0.5235  
}

