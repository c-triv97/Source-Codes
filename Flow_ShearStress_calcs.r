vector.is.empty <- function(x) return(length(x) ==0)

plot.flow <- function(x, title = "Flow"){

    p = ggplot(
        data = x, 
        aes(
            x = Time, 
            y = Filtered.Data
        )
    )+
    geom_line(
        color = "darkgrey"
    )+
    labs(
        x = "Time (s)",
        y = "Flow (ml/min)",
        subtitle = title
    )+
    theme_scientific

    print(p)
}

read.files <- function(folder = getwd(), 
                       names = c(),
                       to.plot = FALSE,
                       height = 0.3, 
                       width  = 28
){

    files = list.files(path = folder,
        pattern = "*filtered.csv",
        full.names = TRUE)

    if (vector.is.empty(files)==TRUE) {
    stop("no processed files in folder")
    }
    
    #print(files)
    
    dat.files = lapply(
        files,
        fread
    )

    dat.files = lapply(
        dat.files, 
        function(x){
            names(x)=make.names(names(x),unique = TRUE)
            return(x)
        }
    )

    if(vector.is.empty(names)==TRUE){
        names = files
    }else{
        names = names
    }

    names(dat.files) = names

    if(to.plot == TRUE){

        lapply(
            names(dat.files), 
            function(x){
                plot.flow(dat.files[[x]],
                          title = x)
    }
    )
    }

    # shear stress

    dat.files = lapply(dat.files, function(x){
        x %>% 
        mutate(
            flow_mm3.s = (Filtered.Data * 1000/60),
            shear = flow_mm3.s/width/(height^2)*6*0.00072, 
            shear.stress = (shear*10))
    })

    return(dat.files)
}

shear_calcs <- function(x,k=30){

    dat <- x %>%
        reframe(
            sample = Sample.Number,
            Time = Time, 
            Instant.Flow = Filtered.Data,
            Instant.Shear.Pa = shear,
            Instant.Shear.Dyn = shear.stress,
            roll.flow = roll::roll_mean(Filtered.Data, width=k),
            roll.shear.pa = roll::roll_mean(shear, width=k),
            roll.shear.dyn = roll::roll_mean(shear.stress, width=k),
            roll.flow.sd = roll_sd(Filtered.Data, width=k),
            roll.shearpa.sd = roll_sd(shear, width=k),
            roll.sheardyn.sd = roll_sd(shear.stress, width=k),
        )
        
    return(dat)
}