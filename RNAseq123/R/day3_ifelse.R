###Day 3
#If/else loops

#conditional statements
#there are a variety of ways we can do this


#one way to do this is with tidyverse
library(tidyverse)
library(ggplot2)

#using tinyverse let's make a column conditionally petal length
#creating a new variable isShort within the function
#we use the = because we are within the function
myiris <- mutate(iris, isShort = )

#let's look at the data
iris
#let's look at the summary
summary(iris)

#looking at the function itself
#elecemntA (logical operator) elementB
ifelse(9>6, "yes", "no")

#make a vector
a <- c(1,3,5,7,9)
b <- c(10, 4, 11, 2, 8)

#conditionally print whether a>b
ifelse(a>b, "greater", "less than")

#creating a new variable isShort within the function
#we use the = because we are within the function
#using the mean value to use 3.758
myiris <- mutate(iris, isShort = ifelse(Petal.Length<3.758, T, F))
glimpse(myiris)

#plot and see
#this will count how many of the instances are short vs long
#face_wrap separates the plot based on the variable after the tild~
ggplot(myiris, aes(x=isShort, fill=Species))+
  geom_bar(stat="count")+
  facet_wrap(~Species)

#conditional statements: use curly brackets
#you can do if without else - but not else without if
#always put the else on the same line as the squiqqly bracket
if(9>6){
  print("Greater!")
} else if (9<6){
  print("Less than!")
} else{
  print("Equal")
}

#FOR LOOP - automate code by iterating through numbers
#for(variable in vector({}
#if variable is i that means index
#
for (i in 1:5) {
  print(i)
  #subset it
  print(a[i])
  print(b[i])
}

#OR
for (i in 1:5) {
  print(i)
  #subset it
  #print(a[i])
  #print(b[i])
  x <- paste(i, a[i], b[i], sep = " ")
}

#make it make more sense
for (i in 1:5) {
  print(i)
  #subset it
  #print(a[i])
  #print(b[i])
  if(a[i]>b[i]) {
    x <- paste("Iteration:", i, "Element A:", a[i], "Element B:" b[i], sep = " ")
  }

}


#Take everything you do and writing it into a function itself
#Writing your own functions
#start simple then will work with larger things
#let's make a function that does what mean does
?mean
#we know it takes our objects and averages them

#make our own
#make up your argument in the brackets
myMean <- function(x, print=T)(
  #calculate the mean
  m<-mean(x)
  if(print){
    stmnt <- paste ("The mean is:", m, sep=" ")
    print(stmnt)
  }
  m
)
#so unless you type out print = F then it's always going to print it out for us

#with necessary arg
myMean(a)


#combine for loop with if else and wrap with a function
#point of a for loop is so you don't have to write out the name of each column
#functions makes it versatile
#loop automates for you

#Make a function that
#let's pseudo code this
# OUTPUT: print the mean of each numeric column in a dataframe
#so how can we do that
#1. Save a temporary variable that is the column name
#2. Check if that column contains numbers
#3. if the column contains numbers, then
  # 4. find the mean of the column (and save it as a variable)
  #5. Print that value
#6. Go to the next iteration (go to the next column name and start again)

#this helps you orient yourself so you can break your column down into smaller and smaller steps



