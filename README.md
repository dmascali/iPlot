# iPlot - Interactive Plot
*iPlot* is a matlab tool for plotting 1-D or 2-D matrices interactively. Data is plotted versus row indexes. 

- `iplot(X)` plots the columns in X (or their spectral amplitude) interactively, one at a time. X can be a matrix or a vector. 
- `iplot(X,Y,...)` plots the columns in X,Y..., one above the other for easy comparison.

## Functionality
*iPlot* functionalities are triggered by pressing keys (see keyboard legend below).  Main functionalities include:
- move forward or backward between the columns of the input matrix ["A" and "D" keys]
- plots the spectral amplitude of the current column [ "F" key ]
- various column ordering (*sequential*, *std +*, *std -*, *random*) ["R" key]

<img align="left"  src="https://user-images.githubusercontent.com/41842946/105090502-c5839680-5a9e-11eb-91cf-35b54e835581.jpg" height="500">
