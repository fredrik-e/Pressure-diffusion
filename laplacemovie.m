function M=laplacemovie(S,iters)

%Creates a movie from the Laplace iterations in S
for i=1:length(S)
    imagesc(S{i}(:,2:size(S{i},2)-1)), colormap jet
    text(0.8,0.9,{['Iteration: ' num2str(iters(i))]},'Units','normalized','Color','Red')
    M(i)=getframe;
end

end