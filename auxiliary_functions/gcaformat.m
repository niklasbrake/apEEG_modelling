function gcaformat(graphicsObject,article,fntsz)
    if(nargin<3)
        fntsz = 7;
    end
    if(nargin<2)
        article = true;
    end
    if(nargin==0)
        nbFormat(gca,article,fntsz)
    elseif(strcmp(graphicsObject.Type,'axes'))
        nbFormat(graphicsObject,article,fntsz)
    elseif(strcmp(graphicsObject.Type,'figure'))
        axs = get(graphicsObject,'children');
        for i = 1:length(axs)
            if(strcmp(axs(i).Type,'axes'))
                nbFormat(axs(i),article,fntsz);
            end
        end
    else
        error('invalid input type');
    end

end
function nbFormat(ax,article,fntsz)
    if(article)
        nbFormat_article(ax,fntsz);
    else
        nbFormat_presentation(ax);
    end
end
function nbFormat_article(ax,fntsz)
    set(ax,'box','off');
    set(ax,'tickdir','out');
    set(ax,'linewidth',0.75);
    set(ax,'fontsize',fntsz);
    set(ax.Title,'FontSize',fntsz);
    xax = get(ax,'xaxis');
    xax.Label.FontSize = fntsz;
    xax.TickLabelRotation = 0;
    yax = get(ax,'yaxis');
    for i = 1:length(yax)
        yax(i).Label.FontSize = fntsz;
        yax(i).TickLabelRotation = 0;
    end
    tickLength = 0.05; % 1/2 mm
    U = ax.Units;
    set(ax,'Units','centimeters');
    L = max(ax.Position(3:4));
    t = tickLength/L;
    set(ax,'TickLength',[t,t]);
    set(ax,'Units',U);
end
function nbFormat_presentation(ax)
    set(ax,'box','off');
    set(ax,'tickdir','out');
    set(ax,'linewidth',1.5);
    set(ax,'fontsize',18);
    set(ax.Title,'FontSize',22);
    xax = get(ax,'xaxis');
    xax.Label.FontSize = 22;
    yax = get(ax,'yaxis');
    for i = 1:length(yax)
        yax(i).Label.FontSize = 22;
    end
    tickLength = 0.2; % 2 mm
    U = ax.Units;
    set(ax,'Units','centimeters');
    L = max(ax.Position(3:4));
    t = tickLength/L;
    set(ax,'TickLength',[t,t]);
    set(ax,'Units',U);
end