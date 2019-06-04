load('locations.RData')

loc_file = Locations[, c('location_code', 'landuse_class', 'block_number', 'xcor', 'ycor')]
loc_file$landuse_class[which(loc_file$landuse_class == '')] = 'HOUSE'

write.table(
  loc_file, file = 'locations.txt',
  quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)



load('repSize.RData')
names(repSize$nb$par) = c('n', 'p')
write.table(
  repSize$nb$par, file = 'repSize.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)


## propLoc
load('propLoc.RData')

# complex
write.table(
  propLoc$L9$par, file = 'propLoc.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)

# simple
types = c(unique(Locations$landuse_class), 'FI')
types = types[-which(types == '')]
props = sapply(sort(types), function(tt)
  if(tt == 'FI'){as.numeric(propLoc$L11$par[tt])}else
  sum(Locations$landuse_class == tt) / (sum(Locations$landuse_class %in% types)) * (1 - as.numeric(propLoc$L11$par['FI'])))
write.table(
  props,
  file = 'propLoc_simple.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)


## where
load('where.RData')

# complex
whereOld = sapply(where$L8, function(tt) tt$par)
whereNew = as.vector(whereOld)
names(whereNew) = sort(c(
  paste(row.names(t(whereOld)), '_a', sep = ''),
  paste(row.names(t(whereOld)), '_b', sep = '')))
write.table(
  whereNew, file = 'where.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)

# simple
types = unique(Locations$landuse_class)
types = types[-which(types == '')]
eval(parse(text = paste('whereNew = 
  c(', paste(sapply(sort(types), function(tt) paste(tt, '_a = 1e-11, ', tt, '_b = 1', sep = '')),
  collapse = ', '), ')', sep = '')))
write.table(
  whereNew, file = 'where_simple.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)


## bins
write.table(
  c(max = 20000, interval = 20), file = 'bins.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)
  

## timeHome
load('timeHome.RData')
timeHomeCoefs = coef(timeHome$I)[1 : 4]
names(timeHomeCoefs) = 
  c('freqMean', 'durnMean', 'freqStd', 'durnStd')
write.table(
  timeHomeCoefs, file = 'timeHome.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)


## freqDurn
load('freqDurn.RData')

# out of town
names(freqDurn$FI$best$par) =
  c('FI_freqMeanA', 'FI_freqStd', 'FI_durnMeanA', 'FI_durnStd')
write.table(
  freqDurn$FI$best$par, file = 'timeFI.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)

# complex
toWrite = numeric()
for(tt in 1 : 10){
  default_values <- c(1, 1, 0, 1, 1, 1, 0, 1, 0)
  names_values <- paste(c("freq.mean.a", "freq.mean.b", "freq.mean.c", "freq.std",
                    "durn.mean.a", "durn.mean.b", "durn.mean.c", "durn.std",
                    "corr"), names(freqDurn$L10)[tt], sep = '.')
  default_values[which(names_values %in% names(freqDurn$L10[[tt]]$best$par))] =
    freqDurn$L10[[tt]]$best$par
  names(default_values) <- paste(names(freqDurn$L10)[tt], c("freqMeanA", "freqMeanB", "freqMeanC", "freqStd",
                    "durnMeanA", "durnMeanB", "durnMeanC", "durnStd",
                    "corr"), sep = '_')
  toWrite = c(toWrite, default_values)
}
write.table(
  toWrite, file = 'freqDurn.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)

# simple
toWrite = numeric()
types = unique(Locations$landuse_class)
types = types[-which(types == '')]
for(tt in sort(types)){
  eval(parse(text = paste('
  values = c(',
    tt, '_freqMeanA = as.numeric(freqDurn$L1[[1]]$NNN$par[1]),',
    tt, '_freqMeanB = 1,',
    tt, '_freqMeanC = 0,',
    tt, '_freqStd = as.numeric(freqDurn$L1[[1]]$NNN$par[2]),',
    tt, '_durnMeanA = as.numeric(freqDurn$L1[[1]]$NNN$par[3]),',
    tt, '_durnMeanB = 1,',
    tt, '_durnMeanC = 0,',
    tt, '_durnStd = as.numeric(freqDurn$L1[[1]]$NNN$par[4]),',
    tt, '_corr = 0)', sep = '')))
  toWrite = c(toWrite, values)
}
write.table(
  toWrite, file = 'freqDurn_simple.txt',
  quote = FALSE, sep = '\t', row.names = TRUE, col.names = FALSE)




freqDurn$L1[[1]]$NNN$par[1]


prop = read.csv('kellyData.csv')
prop<-subset(prop, interview_number<=8)
prop<-subset(prop, entrances != "NULL")
prop<-subset(prop, entrances!= "NA")
prop<-subset(prop, hours_in_house != "NULL")
prop<-subset(prop, hours_in_house != "NA")
prop<-subset(prop, dtr != "NULL")
prop<-subset(prop, dtr != "NA")

prop$loc_int = sapply(1 : nrow(prop), function(rr)
  paste(prop[rr, c('location_code', 'interview_number')], collapse = '_'))


write.table('', file = 'houseComp.txt',
  append = FALSE, quote = FALSE, eol = '', row.names = FALSE, col.names = FALSE)
for(ii in unique(prop$loc_int)){
  write.table(
    cbind(t(subset(prop, loc_int == ii)[, c('hours_in_house', 'body_surface')]), c(-1, -1)),
    file = 'houseComp.txt',
    append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}




prop$census = paste(prop$location_code,prop$interview_number,sep='_')

censuses = data.frame(
  census = names(table(prop$census)),
  numres = as.numeric(table(prop$census)))
censuses = censuses[-which(censuses$census=='MY253A_4'),]

locs.sim = read.csv('locations_20150801.csv')
p = optimize(f=function(p){
  abs(sum(dgeom(censuses$numres,p)/sum(dgeom(censuses$numres,p))*censuses$numres)*table(locs.sim$landuse)['HOUSE']-200000)},
  interval=c(0,1))$minimum
syn.locs = sample(nrow(censuses),table(locs.sim$landuse)['HOUSE'],replace=T,prob=dgeom(censuses$numres,p))

census.real = read.csv('census.csv')
total.syn.pops = sum(censuses$numres[syn.locs])
males.real = census.real$Male / sum(census.real$Total.1) * total.syn.pops
females.real = census.real$Female / sum(census.real$Total.1) * total.syn.pops
males.tally = ceiling(males.real)
females.tally = ceiling(females.real)

write.table('', file = 'wholePop.txt',
            append = FALSE, quote = FALSE, eol = '', row.names = FALSE, col.names = FALSE)
for(ii in syn.locs){
  sex.write = prop$sex[prop$census==censuses$census[ii]]
  age.write = prop$age_at_interview[prop$census==censuses$census[ii]]
  for(jj in 1:censuses$numres[ii]){
    if(sex.write[jj] == 'F' & sum(females.tally) > 0){
      if(females.tally[age.write[jj]+1] > 0){
        females.tally[age.write[jj]+1] = females.tally[age.write[jj]+1] - 1
      }else{
        if(age.write[jj] < 18 & sum(females.tally[1:18]) > 0){
          age.write[jj] = sample(0:17,1,prob=females.tally[1:18])
          females.tally[age.write[jj]+1] = females.tally[age.write[jj]+1] - 1
        }else{
          age.write[jj] = sample(18:98,1,prob=females.tally[19:99])
          females.tally[age.write[jj]+1] = females.tally[age.write[jj]+1] - 1          
        }
      }
    }else{
      sex.write[jj] = 'M'
      if(males.tally[age.write[jj]+1] > 0){
        males.tally[age.write[jj]+1] = males.tally[age.write[jj]+1] - 1
      }else{
        if(age.write[jj] < 18 & sum(males.tally[1:18]) > 0){
          age.write[jj] = sample(0:17,1,prob=males.tally[1:18])
          males.tally[age.write[jj]+1] = males.tally[age.write[jj]+1] - 1
        }else{
          age.write[jj] = sample(18:98,1,prob=males.tally[19:99])
          males.tally[age.write[jj]+1] = males.tally[age.write[jj]+1] - 1          
        }
      }
    }
  }
  
  write.table(
    paste(
      paste(
        sex.write,
        age.write,
        sep='\t',
        collapse='\t'),
      'N', sep='\t'),
    file = 'wholePop.txt',
    append = TRUE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}
