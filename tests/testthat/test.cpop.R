

context("testing cpop")

test_that("test that cpop predicts the correct RSS and changepoints using default parameter values",
{
   load(file="test.1.RData")
   cpop.res<-cpop(y,x)
   cpop.RSS<-sum(fitted(cpop.res)$RSS)
   expect_equal(cpop.RSS,RSS)
   expect_equal(changepoints(cpop.res)$location,x[out$changepoints[2:4]])
})

test_that("test that cpop predicts the correct RSS and changepoints using non default beta value",
{
   load(file="test.2.RData")
   cpop.res<-cpop(y,x,beta=2)
   cpop.RSS<-sum(fitted(cpop.res)$RSS)
   expect_equal(cpop.RSS,RSS)
   expect_equal(changepoints(cpop.res)$location,x[out$changepoints[2:24]])
})





test_that("test that fitted predicts the correct RSS",
{
   load(file="test.3.RData")
   cpop.res<-cpop(y,x,sd=1)   
   expect_equal(sum(fitted(cpop.res)$RSS),RSS)
})


test_that("test that fitted predicts the correct RSS",
{
   load(file="test.4.RData")
   cpop.res<-cpop(y,x,sd=1)
   expect_equal(sum(fitted(cpop.res)$RSS),RSS)
})


test_that("test that results from fitted can be used to calculate the cost correctly",
{
   load(file="test.5.RData")
   cpop.res<-cpop(y,x,sd=1)
   expect_equal(cost(cpop.res),cost)
})



test_that("test default value of sd",
{
   load(file="test.6.RData")
   cpop.res<-cpop(y,x)
   expect_equal(cost(cpop.res),cost)
})


test_that("test non default values of sd",
{
   load(file="test.7.RData")
   cpop.res<-cpop(y,x,sd=sqrt(2))
   expect_equal(cost(cpop.res),cost)
   load(file="test.8.RData")
   cpop.res<-cpop(y,x,sd=2)
   expect_equal(cost(cpop.res),cost)
})


test_that("test for the effects of setting minseglen greater than shortest distance between changepoints",
{
   load(file="test.9.RData")
   cpop.res<-cpop(y,x,sd=2)
   expect_equal(cost(cpop.res),cost)
   cpop.minseglen.res<-cpop(y,x,sd=2,minseglen=22)
   expect_equal(cost(cpop.res),cost(cpop.minseglen.res))
   cpop.minseglen.res<-cpop(y,x,sd=2,minseglen=23)
   expect_false(isTRUE(all.equal(cost(cpop.res),cost(cpop.minseglen.res))))
   cpop.minseglen.res<-cpop(y,x,sd=2,minseglen=26)
   expect_equal(changepoints(cpop.minseglen.res)$location[1],110)
})


test_that("test use of non default value for grid",
{
   set.seed(1)
   x<-1:200
   changepoints<-c(0,25,50,100)
   change.slope<-c(0.2,-0.3,0.2,-0.1)
   y<-simulate(x,changepoints,change.slope,1)
   cpop.res<-cpop(y,x,sd=2)
   cpop.grid.res<-cpop(y,x,grid=c(0.5,1.5,99.0),sd=2)
   expect_false(isTRUE(all.equal(cost(cpop.res),cost(cpop.grid.res))))
   expect_equal(changepoints(cpop.grid.res)$location[1],99)
})



test_that("test use of non unit locations data",
{
   load(file="test.10.RData")
   cpop.res <- cpop(y,x,beta=2*log(length(y)),sd=0.1)
   expect_equal(sum(fitted(cpop.res)$RSS),RSS)
})



test_that("test that simulate produces data consistently",
{
   set.seed(1)
   changepoints=c(0,25,50,100)
   change.slope=c(0.2,-0.3,0.2,-0.1)
   x=1:200
   sig=1+x/200
   y<-simulate(x,changepoints,change.slope,sig)
   load("test.11.RData")
   expect_equal(y,test.simulate.y)
})








