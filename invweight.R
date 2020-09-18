## ------------------------
## invweight function
## ------------------------

## Translation from Matlab -> R of invweight function from Keller & Schoene (2012, Nature)
## Geospatial weighting function (analagous to statistical "declustering' methods)
## input latitude, longitude and age data
## output weighting coefficient "k" for geospatial and temporal weighting 

invweight <- function(lat, lon, age){

age_bin = 38  
k <- numeric(length(lat))
p=2
for(i in 1:length(lat)){{
    if(is.na(lat) | is.na(lon) | is.na(age)){
    k[i]=Inf  
    }else{
    k[i]=sum((1/((180/pi*acos(sin(lat[i]*pi/180)*sin(lat*pi/180)+cos(lat[i]*pi/180)*cos(lat*pi/180)*cos(lon[i]*pi/180-lon*pi/180))/1.8)^p+1)+1/((((age[i])-age)/age_bin)^p+1)), na.rm=TRUE) 
  
    ## Matlab code from Keller & Schoene (2012, Nature)
    #___________
    ## k(i)=nansum(1./((180/pi*acos(sin(lat(i)*pi/180).*sin(lat*pi/180)+cos(lat(i)*pi/180).*cos(lat*pi/180).*cos(lon(i)*pi/180-lon*pi/180))/1.8).^p+1)...
    ## +1./((((age(i))-age)/38).^p+1));        
    
    
if(k[i] <0){ # Check for anomalous k values
  print(i)
  print(lat[i])
  print(lon[i])
  print(age[i])
}
}      
}}
return(k)
}
 
