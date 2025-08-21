Created by shazia
first mini-task of bioinformatics 1 in 4th semester which is now slightly updated by me

Reads two input strings and sets scoring (match +2, mismatch −1, gap −1).

Fills a DP score matrix H: for each cell it takes max(0, diag+score, up+gap, left+gap) and remembers the direction.

Starts traceback from the highest‐scoring cell to rebuild the best local alignment (adding gaps where needed).

Prints the score matrix, the best local score, and the two aligned substrings with a match guide
