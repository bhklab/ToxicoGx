### TODO:: Add updating of sensitivity Number tables
#' A function to update drug ids
#' @examples
#' data(TGGATESsmall)
#' updateDrugId(TGGATESsmall, new.ids = drugNames(TGGATESsmall))
#'
#' @param tSet [object] A ToxicoSet object to be updates
#' @param new.ids [character] A character vector of ids to update with
#'
#' @return \code{none} Updates the drug ids in the ToxicoSet
#'
#' @keywords internal
#' @export
updateDrugId <- function(tSet, new.ids = vector("character")){

    if (length(new.ids)!= nrow(drugInfo(tSet))){
        stop("Wrong number of drug identifiers")
    }

    if(datasetType(tSet)=="sensitivity" | datasetType(tSet)=="both"){
        myx <- match(sensitivityInfo(tSet)[,"drugid"],rownames(drugInfo(tSet)))
        sensitivityInfo(tSet)[,"drugid"] <- new.ids[myx]

    }
    if(datasetType(tSet)=="perturbation"|datasetType(tSet)=="both"){
        molecularProfilesSlot(tSet) <- lapply(molecularProfilesSlot(tSet), function(SE){

            myx <- match(SummarizedExperiment::colData(SE)[["drugid"]],rownames(drugInfo(tSet)))
            SummarizedExperiment::colData(SE)[["drugid"]]  <- new.ids[myx]
            return(SE)
        })
    }


    if(any(duplicated(new.ids))){
        warning("Duplicated ids passed to updateDrugId. Merging old ids into the
            same identifier")

        if(ncol(sensNumber(tSet)) > 0){
            sensMatch <- match(colnames(sensNumber(tSet)), rownames(drugInfo(tSet)))
        }
        if(dim(pertNumber(tSet))[[2]] > 0){
            pertMatch <- match(dimnames(pertNumber(tSet))[[2]], rownames(drugInfo(tSet)))
        }
        curMatch <- match(rownames(tSet@curation$drug),rownames(drugInfo(tSet)))

        duplId <- unique(new.ids[duplicated(new.ids)])
        for(id in duplId){

            if (ncol(sensNumber(tSet))>0){
                myx <- which(new.ids[sensMatch] == id)
                sensNumber(tSet)[,myx[1]] <- apply(sensNumber(tSet)[,myx], 1, sum)
                sensNumber(tSet) <- sensNumber(tSet)[,-myx[-1]]
                # sensMatch <- sensMatch[-myx[-1]]
            }
            if (dim(pertNumber(tSet))[[2]]>0){
                myx <- which(new.ids[pertMatch] == id)
                pertNumber(tSet)[,myx[1],] <- apply(pertNumber(tSet)[,myx,], c(1,3), sum)
                pertNumber(tSet) <- pertNumber(tSet)[,-myx[-1],]
                # pertMatch <- pertMatch[-myx[-1]]
            }

            myx <- which(new.ids[curMatch] == id)
            tSet@curation$drug[myx[1],] <- apply(tSet@curation$drug[myx,], 2, paste, collapse="///")
            tSet@curation$drug <- tSet@curation$drug[-myx[-1],]
            # curMatch <- curMatch[-myx[-1]]

            myx <- which(new.ids == id)
            drugInfo(tSet)[myx[1],] <- apply(drugInfo(tSet)[myx,], 2, paste, collapse="///")
            drugInfo(tSet) <- drugInfo(tSet)[-myx[-1],]
            new.ids <- new.ids[-myx[-1]]
            if(ncol(sensNumber(tSet))>0){
                sensMatch <- match(colnames(sensNumber(tSet)), rownames(drugInfo(tSet)))
            }
            if(dim(pertNumber(tSet))[[2]]>0){
                pertMatch <- match(dimnames(pertNumber(tSet))[[2]], rownames(drugInfo(tSet)))
            }
            curMatch <- match(rownames(tSet@curation$drug),rownames(drugInfo(tSet)))
        }
    } else {
        if (dim(pertNumber(tSet))[[2]]>0){
            pertMatch <- match(dimnames(pertNumber(tSet))[[2]], rownames(drugInfo(tSet)))
        }
        if (ncol(sensNumber(tSet))>0){
            sensMatch <- match(colnames(sensNumber(tSet)), rownames(drugInfo(tSet)))
        }
        curMatch <- match(rownames(tSet@curation$drug),rownames(drugInfo(tSet)))
    }

    if (dim(pertNumber(tSet))[[2]]>0){
        dimnames(pertNumber(tSet))[[2]] <- new.ids[pertMatch]
    }
    if (ncol(sensNumber(tSet))>0){
        colnames(sensNumber(tSet)) <- new.ids[sensMatch]
    }
    rownames(tSet@curation$drug) <- new.ids[curMatch]
    rownames(drugInfo(tSet)) <- new.ids


    return(tSet)
}