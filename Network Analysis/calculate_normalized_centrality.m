function nc = calculate_normalized_centrality(G, cName)

     c = centrality(G, cName);
     N = numnodes(G);
     
     if strcmp(cName, 'closeness') || strcmp(cName, 'Closeness')
         nc = c * (N - 1);
     elseif strcmp(cName, 'betweenness') || strcmp(cName, 'Betweenness')
         nc = 2*c / ( (N-1)*(N-2) );
     else
         nc = c;
     end

end