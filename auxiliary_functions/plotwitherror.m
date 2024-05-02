function varargout =plotwitherror(x,y,isstdrr,varargin)

	if(strcmp(isstdrr,'CI') | isstdrr == true)
		ymu = nanmean(y,2);
		z = icdf('norm',0.975,0,1);
		ylo = ymu-z*stderror(y')';
		yhi = ymu+z*stderror(y')';
	elseif(strcmp(isstdrr,'SE') | isstdrr == false)
		ymu = nanmean(y,2);
		ylo = ymu-std(y')';
		yhi = ymu+std(y')';
	elseif(strcmp(isstdrr,'Q'))
		ymu = nanmedian(y,2);
		ylo = quantile(y,0.025,2);
		yhi = quantile(y,0.925,2);
	elseif(strcmp(isstdrr,'M'))
		ymu = nanmean(y,2);
		ylo = min(y,[],2);
		yhi = max(y,[],2);
	end

	h = plot(x,ymu,varargin{:}); hold on;

	idcs = ~isnan(ymu);
	iRegion = 1;
	changeRegions = [0;(diff(idcs)==-1)];
	region=cell(1);
	for i = 1:length(idcs)
		if(idcs(i)==1)
			region{iRegion} = [region{iRegion},i];
		else
			if(changeRegions(i))
				iRegion = iRegion+1;
				region{iRegion} = [];
			end
		end
	end

	for i = 1:length(region)
		xfvec = [x(region{i}),flip(x(region{i}))];
		yfvec = [ylo(region{i});flip(yhi(region{i}))];
		F = fill(xfvec(:),yfvec(:),'k');

		if(isempty(F))
			continue;
		end
		F.FaceAlpha = 0.2;
		F.FaceColor = h.Color;
		F.EdgeColor = 'none';
	end

	if(nargout>0)
		varargout = {h};
	end
end
function SE = stderror(X)
	SE = nanstd(X,0,1)./sqrt(sum(~isnan(X),1));
end