# Load Synapse CGH workflow
loadWorkflow <- function(){
  cat('--------------------------\n')
  cat('-  Loading CGH Workflow  -\n')
  cat('--------------------------\n')
  # Required R packages
  code <- synGet('syn2121335'); source(code@filePath)
  cat('Class definitions...\t')
  code <- synGet('syn2121329'); source(code@filePath)
  cat('Done.\nGenerics definitions...\t')
  code <- synGet('syn2121330'); source(code@filePath)
  cat('Done.\nAccessor functions...\t')
  code <- synGet('syn2121328'); source(code@filePath)
  cat('Done.\nShow methods...\t')
  code <- synGet('syn2121336'); source(code@filePath)
  cat('Done.\nMethods...\t')
  code <- synGet('syn2121332'); source(code@filePath)
  cat('Done.\nHelper functions...\n')
  code <- synGet('syn2121331'); source(code@filePath)
  cat('Done.\nSynapse specific functions, including shiny builders...\t')
  code <- synGet('syn2121337'); source(code@filePath)
  cat('Done.\n')
}
loadWorkflow()