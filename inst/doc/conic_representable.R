## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6,
  eval = TRUE,
  error = TRUE
)

## ----conic, echo = FALSE, fig.cap = cap_1, figure = TRUE----------------------
cap_1 <- "Simple decision tree to evaluate if a given optimization problem can be solved with a conic optimization solver."
library("DiagrammeR")
mermaid('
%%{init: { "theme": "neutral" } }%%
graph TB
    A["Can my optimization problem be solved via conic optimization?"] --> B{"Convex?"};
    B -- yes --> C{"Representable?"};
    B -- no --> D["Not solvable (Exceptions A)."];
    C -- yes --> E["Solvable via conic optimization."];
    C -- no --> F["Not solvable (Exceptions B)."];
    click D "#ExceptionsA" "For more information see Section `Exceptions A`"
    click F "#ExceptionsB" "For more information see Section `Exceptions B`"
')

