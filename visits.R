### Code from http://www.tatvic.com/blog/calender-heatmap-with-google-analytics-data/

#Load RGoogleAnalytics library
library("RGoogleAnalytics")

# Create query builder object
query <- QueryBuilder()

# Authorize your account and paste the accesstoken
access_token <- query$authorize()

# Create a new Google Analytics API object
ga <- RGoogleAnalytics()
ga.profiles <- ga$GetProfileData(access_token)

# List the GA profiles
ga.profiles  # select index corresponds to your profile and set it to query string

# For example if index is 7 of your GA profile then set ga.profiles$id[7] in
# query$Init() method given below
# Build the query string
query$Init(start.date = "2012-12-01", # Set start date
           end.date = "2013-03-17", # Set end date
           dimensions = "ga:date",
           metrics = "ga:visits,ga:transactions",
           max.results = 10000,
           table.id = paste("ga:",ga.profiles$id[1],sep="",collapse=","),
           access_token=access_token)

# Make a request to get the data from the API
ga.data <- ga$GetReportData(query)  # data will be stored in this data frame

# Set date in format YYYY-MM-DD (to use into heatmap calender)
ga.data$date <- as.Date(as.character(ga.data$date),format="%Y%m%d")


# Recommended R version - 2.15.1 or higher
# install required  library by using the command install.packages('libraryname”)
# For example install.packages('ggplot2”)
# Required library
library('quantmod')
library('ggplot2')
library('reshape2')
library('plyr')
library('scales')

# Set extracted data  to this data frame
data <-  ga.data

# Run commands listed below
data$year <- as.numeric(as.POSIXlt(data$date)$year+1900)
data$month <- as.numeric(as.POSIXlt(data$date)$mon+1)
data$monthf <- factor(data$month,levels=as.character(1:12),
                      labels=c("Jan","Feb","Mar","Apr","May","Jun",
                               "Jul","Aug","Sep","Oct","Nov","Dec"),
                      ordered=TRUE)
data$weekday <- as.POSIXlt(data$date)$wday
data$weekdayf <- factor(data$weekday,levels=rev(0:6),
                        labels=rev(c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")),
                        ordered=TRUE)
data$yearmonth <- as.yearmon(data$date)
data$yearmonthf <- factor(data$yearmonth)
data$week <- as.numeric(format(as.Date(data$date),"%W"))
data <- ddply(data,.(yearmonthf),transform,monthweek=1+week-min(week))

# Plot for visits
P_visits <- ggplot(data, aes(monthweek, weekdayf, fill = visits)) +
  geom_tile(colour = "white") +
  facet_grid(year~monthf) +
  scale_fill_gradient(high="#D61818",low="#B5E384") +
  labs(title = "Time-Series Calendar Heatmap") +
  xlab("Week of Month") +
  ylab("")

# View plot
P_visits