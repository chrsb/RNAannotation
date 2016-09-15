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

# getting only the first 5 k entries
data <- data[1:2500, ]
classes <- classes[1:2500, ]
# it's doing something funky if the name column is left in (some sort of giant data structure gets built and I end up with 'Error: cannot allocate vector of size 46.0 Gb')
data$V1 <- NULL

# glueing the class column to the data
data <- data.frame(data, Class = classes[, 2])

set.seed(100)

trainingClasses <- createDataPartition(data$Class, p = .75, list = FALSE)
trainingSet <- data[trainingClasses, ]
testSet <- data[-trainingClasses, ]
fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10)
						   
set.seed(100)

results <- 1:nrow(config)

for (i in 1:nrow(config)) {
    # checks if the second config column is empty and parses the grid function if it's there
    if (as.character(config[i, 2] != "")) {
        #results[i] <- train(Class ~ ., data = trainingSet, 
        result <- train(Class ~ ., data = trainingSet, 
                          method = as.character(config[i, 1]), 
                          trControl = fitControl,
                          ## This last option is actually one
                          ## for gbm() that passes through
                          verbose = FALSE,
                          tuneGrid = eval(parse(text=as.character(config[i, 2]))))
        print(result) 
        plot(result) 
        
        # this needs to move to its own loop when I figure out how to accumulate results in a list
        x_test <- testSet[, 1:7]
        y_test <- testSet[, 8]
        predictions <- predict(result, x_test)
        confusionMatrix(predictions, y_test)
    } else {
        #results[i] <- train(Class ~ ., data = trainingSet, 
        result <- train(Class ~ ., data = trainingSet, 
                          method = as.character(config[i, 1]), 
                          trControl = fitControl,
                          ## This last option is actually one
                          ## for gbm() that passes through
                          verbose = FALSE)
        print(result) 
        plot(result) 
        
        # this needs to move to its own loop when I figure out how to accumulate results in a list
        x_test <- testSet[, 1:7]
        y_test <- testSet[, 8]
        predictions <- predict(result, x_test)
        confusionMatrix(predictions, y_test)
    }
}

for (i in 1:nrow(config)) {
    # predict - the prediction stuff moves here
}

