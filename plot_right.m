function [ax2,hpl]=plot_right(hax,xval,yval,varargin)
%[ax2,hpl]=plot_right(hax,xval,yval)
%[ax2,hpl]=plot_right(hax,xval,yval,linespec)

ax2 = axes('Position',get(hax,'Position'),...
        'XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','k',...
        'XMinorGrid','off','YGrid','off');
    hold on
    hpl=plot(xval,yval,varargin{:},'Parent',ax2);
    
    set(ax2,'Position',get(hax,'Position'),...
        'XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor','k','XLim',get(hax,'XLim'));
    linkaxes([hax,ax2],'x');
    ch=get(get(hax,'Parent'),'Children');
    rind= hax==ch;
    pfind= ax2==ch;
    ch(rind)=ax2;
    ch(pfind)=hax(1);
    %set(get(hax,'Parent'),'Children',ch)