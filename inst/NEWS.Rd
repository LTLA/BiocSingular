\name{BiocSingularnews}
\title{BiocSingular News}
\encoding{UTF-8}

\section{Version 1.2.0}{\itemize{
\item Added the ResidualMatrix class for computing PCA on residuals efficiently.

\item Fixed runIrlba() to avoid errors at the limit of available PCs.

\item Added the FastAutoParam class to automatically choose a fast SVD depending on the matrix representation.
}}

\section{Version 1.0.0}{\itemize{
\item
New package BiocSingular, for performing SVD and PCA within a Bioconductor framework.
}}
