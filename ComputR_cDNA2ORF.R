fetch_endpoint <- function(server, request, content_type){
  
  r <- GET(paste(server, request, sep = ""), accept(content_type))
  
  stop_for_status(r)
  
  if (content_type == 'application/json'){
    return (fromJSON(content(r, "text", encoding = "UTF-8")))
  } else {
    return (content(r, "text", encoding = "UTF-8"))
  }
}

fetch_endpoint_POST <- function(server, request, data, content_type='application/json'){
  
  r <- POST(paste(server, request, sep = ""), content_type("application/json"), accept("application/json"), body = data)
  
  stop_for_status(r)
  
  if (content_type == 'application/json'){
    return (fromJSON(content(r, "text", encoding = "UTF-8")))
  } else {
    return (content(r, "text", encoding = "UTF-8"))
  }
}

findPotentialStartsAndStops <- function(sequence)
{ codons <- c("ATG", "TAA", "TAG", "TGA")
# Find the number of occurrences of each type of potential start or stop codon
for (i in 1:4){
  codon <- codons[i]
  # Find all occurrences of codon "codon" in sequence "sequence"
  occurrences <- matchPattern(codon, sequence)
  # Find the start positions of all occurrences of "codon" in sequence "sequence"
  codonpositions <- attr(occurrences@ranges,"start")
  # Find the total number of potential start and stop codons in sequence "sequence"
  numoccurrences <- length(codonpositions)
  if (i == 1)
  {
    # Make a copy of vector "codonpositions" called "positions"
    positions <- codonpositions
    # Make a vector "types" containing "numoccurrences" copies of "codon"
    types <- rep(codon, numoccurrences)
  }
  else
  {
    # Add the vector "codonpositions" to the end of vector "positions":
    positions <- append(positions, codonpositions, after=length(positions))
    # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
    types <- append(types, rep(codon, numoccurrences), after=length(types))
  }
}
# Sort the vectors "positions" and "types" in order of position along the input sequence:
indices <- order(positions)
positions <- positions[indices]
types <- types[indices]
# Return a list variable including vectors "positions" and "types":
mylist <- list(positions,types)
return(mylist)
}

plotPotentialStartsAndStops <- function(sequence)
{
  # Define a vector with the sequences of potential start and stop codons
  codons <- c("ATG", "TAA", "TAG", "TGA")
  # Find the number of occurrences of each type of potential start or stop codon
  for (i in 1:4)
  {
    codon <- codons[i]
    # Find all occurrences of codon "codon" in sequence "sequence"
    occurrences <- matchPattern(codon, sequence)
    # Find the start positions of all occurrences of "codon" in sequence "sequence"
    codonpositions <- attr(occurrences@ranges,"start")
    # Find the total number of potential start and stop codons in sequence "sequence"
    numoccurrences <- length(codonpositions)
    if (i == 1)
    {
      # Make a copy of vector "codonpositions" called "positions"
      positions   <- codonpositions
      # Make a vector "types" containing "numoccurrences" copies of "codon"
      types       <- rep(codon, numoccurrences)
    }
    else
    {
      # Add the vector "codonpositions" to the end of vector "positions":
      positions   <- append(positions, codonpositions, after=length(positions))
      # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
      types       <- append(types, rep(codon, numoccurrences), after=length(types))
    }
  }
  # Sort the vectors "positions" and "types" in order of position along the input sequence:
  indices <- order(positions)
  positions <- positions[indices]
  types <- types[indices]
  # Make a plot showing the positions of the start and stop codons in the input sequence:
  # Draw a line at y=0 from 1 to the length of the sequence:
  x  <- c(1,nchar(sequence))
  y <- c(0,0)
  plot <- plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide", ylab="Reading frame",
       main="Predicted start (red) and stop (blue) codons")
  segments(1,1,nchar(sequence),1)
  segments(1,2,nchar(sequence),2)
  # Add the x-axis at y=0:
  axis(1, pos=0)
  # Add the y-axis labels:
  text(0.9,0.5,"+1")
  text(0.9,1.5,"+2")
  text(0.9,2.5,"+3")
  # Draw in each predicted start/stop codon:
  numcodons <- length(positions)
  for (i in 1:numcodons)
  {
    position <- positions[i]
    type <- types[i]
    remainder <- (position-1) %% 3
    if    (remainder == 0) # +1 reading frame
    {
      if (type == "ATG") { segments(position,0,position,1,lwd=1,col="red") }
      else               { segments(position,0,position,1,lwd=1,col="blue")}
    }
    else if (remainder == 1)
    {
      if (type == "ATG") { segments(position,1,position,2,lwd=1,col="red") }
      else               { segments(position,1,position,2,lwd=1,col="blue")}
    }
    else if (remainder == 2)
    {
      if (type == "ATG") { segments(position,2,position,3,lwd=1,col="red") }
      else               { segments(position,2,position,3,lwd=1,col="blue")}
    }
  }
  return(plot)
}

findORFsinSeq <- function(sequence)
{
  require(Biostrings)
  # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
  mylist <- findPotentialStartsAndStops(sequence)
  positions <- mylist[[1]]
  types <- mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts <- numeric()
  orfstops <- numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths <- numeric()
  # Print out the positions of ORFs in the sequence:
  # Find the length of vector "positions"
  numpositions <- length(positions)
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2)
  {
    for (i in 1:(numpositions-1))
    {
      posi <- positions[i]
      typei <- types[i]
      found <- 0
      while (found == 0)
      {
        for (j in (i+1):numpositions)
        {
          posj  <- positions[j]
          typej <- types[j]
          posdiff <- posj - posi
          posdiffmod3 <- posdiff %% 3
          # Add in the length of the stop codon
          orflength <- posj - posi + 3
          if (typei == "ATG" && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0)
          {
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs <- length(orfstops)
            usedstop <- -1
            if (numorfs > 0)
            {
              for (k in 1:numorfs)
              {
                orfstopk <- orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop <- 1 }
              }
            }
            if (usedstop == -1)
            {
              orfstarts <- append(orfstarts, posi, after=length(orfstarts))
              orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths <- append(orflengths, orflength, after=length(orflengths))
            }
            found <- 1
            break
          }
          if (j == numpositions) { found <- 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices <- order(orfstarts)
  orfstarts <- orfstarts[indices]
  orfstops <- orfstops[indices]
  # Find the lengths of the ORFs that we have
  orflengths <- numeric()
  numorfs <- length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart <- orfstarts[i]
    orfstop <- orfstops[i]
    orflength <- orfstop - orfstart + 1
    orflengths <- append(orflengths,orflength,after=length(orflengths))
  }
  mylist <- list(orfstarts, orfstops, orflengths)
  return(mylist)
}

plotORFsinSeq <- function(sequence)
{
  # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
  mylist <- findPotentialStartsAndStops(sequence)
  positions <- mylist[[1]]
  types <- mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts <- numeric()
  orfstops <- numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths <- numeric()
  # Print out the positions of ORFs in the sequence:
  numpositions <- length(positions) # Find the length of vector "positions"
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2)
  {
    for (i in 1:(numpositions-1))
    {
      posi <- positions[i]
      typei <- types[i]
      found <- 0
      while (found == 0)
      {
        for (j in (i+1):numpositions)
        {
          posj <- positions[j]
          typej <- types[j]
          posdiff <- posj - posi
          posdiffmod3 <- posdiff %% 3
          orflength <- posj - posi + 3 # Add in the length of the stop codon
          if (typei == "ATG" && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0)
          {
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs <- length(orfstops)
            usedstop <- -1
            if (numorfs > 0)
            {
              for (k in 1:numorfs)
              {
                orfstopk <- orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop <- 1 }
              }
            }
            if (usedstop == -1)
            {
              orfstarts <- append(orfstarts, posi, after=length(orfstarts))
              orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths <- append(orflengths, orflength, after=length(orflengths))
            }
            found <- 1
            break
          }
          if (j == numpositions) { found <- 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices <- order(orfstarts)
  orfstarts <- orfstarts[indices]
  orfstops <- orfstops[indices]
  # Make a plot showing the positions of ORFs in the input sequence:
  # Draw a line at y=0 from 1 to the length of the sequence:
  x <- c(1,nchar(sequence))
  y <- c(0,0)
  plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide", ylab="Reading frame", main="Predicted ORFs")
  segments(1,1,nchar(sequence),1)
  segments(1,2,nchar(sequence),2)
  # Add the x-axis at y=0:
  axis(1, pos=0)
  # Add the y-axis labels:
  text(0.9,0.5,"+1")
  text(0.9,1.5,"+2")
  text(0.9,2.5,"+3")
  # Make a plot of the ORFs in the sequence:
  numorfs <- length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart <- orfstarts[i]
    orfstop <- orfstops[i]
    remainder <- (orfstart-1) %% 3
    if    (remainder == 0) # +1 reading frame
    {
      rect(orfstart,0,orfstop,1,col="cyan",border="black")
    }
    else if (remainder == 1)
    {
      rect(orfstart,1,orfstop,2,col="cyan",border="black")
    }
    else if (remainder == 2)
    {
      rect(orfstart,2,orfstop,3,col="cyan",border="black")
    }
  }
}

printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)
{
  require(Biostrings)           # This function requires the Biostrings package
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
  starts  <- seq(1, alnlen, by=chunksize)
  n       <- length(starts)
  seq1alnresidues <- 0
  seq2alnresidues <- 0
  for (i in 1:n) {
    chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
    chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
    # Find out how many gaps there are in chunkseq1aln:
    gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
    # Find out how many gaps there are in chunkseq2aln:
    gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
    # Calculate how many residues of the first sequence we have printed so far in the alignment:
    seq1alnresidues <- seq1alnresidues + chunksize - gaps1
    # Calculate how many residues of the second sequence we have printed so far in the alignment:
    seq2alnresidues <- seq2alnresidues + chunksize - gaps2
    if (returnlist == 'FALSE')
    {
      print(paste(chunkseq1aln,seq1alnresidues))
      print(paste(chunkseq2aln,seq2alnresidues))
      print(paste(' '))
    }
  }
  if (returnlist == 'TRUE')
  {
    vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
    vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
    mylist <- list(vector1, vector2)
    return(mylist)
  }
}
