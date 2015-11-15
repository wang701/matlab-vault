function varargout = zpgui(varargin)
% ZPGUI  Zero Pole dragging Graphic User Interface
% Allows you to add, remove, and move the zeros and 
% poles of a filter interactively.
% To begin, type:
%    zpgui

%   Author: Tom Krauss, 9/1/98

global zh ph Lresp Nfft b a z h Y
global ax1 ax2

if nargin == 0
    action = 'init';
else
    action = varargin{1};
end

switch action

case 'init'
    set(0,'defaultaxesfontsize',20)

    subplot(2,2,1)
    z = .4;
    p = .5;
    [zh,ph,cruff]=zplane(z,p);
    %set(cruff,'hittest','off')
    ax1 = gca;
    ylabel('Imaginary Part')
    xlabel('Real Part')
    
    uicontrol('style','pushbutton',...
         'string','Add Zero',...
         'fontsize',14,...
         'units','normalized',...
         'position',[.55 .6 .2 .1],...
         'callback','zpgui(''addzero'')');
    uicontrol('style','pushbutton',...
         'string','Remove Zero',...
         'fontsize',14,...
         'units','normalized',...
         'position',[.8 .6 .2 .1],...
         'callback','zpgui(''removezero'')');
    uicontrol('style','pushbutton',...
         'units','normalized',...
         'fontsize',14,...
         'position',[.55 .8 .2 .1],'string','Add Pole',...
         'callback','zpgui(''addpole'')');
    uicontrol('style','pushbutton',...
         'units','normalized',...
         'fontsize',14,...
         'position',[.8 .8 .2 .1],'string','Remove Pole',...
         'callback','zpgui(''removepole'')');
         
    subplot(2,1,2)

    [b,a]=zp2tf(z,p,1);
    
    Nfft = 128;
    
    Y = fft(b,Nfft)./fft(a,Nfft);
    Y = Y/max(abs(Y));
    Lresp = plot((0:Nfft-1)/Nfft - .5, 20*log10(fftshift(abs(Y))),'linewidth',5);
    ax2 = gca;
    set(ax2,'xlim',[-.5 .5])
    set(ax2,'ylim',[-50 20])
    grid on
    xlabel('Frequency')
    ylabel('Magnitude (dB)')
    
    set(Lresp,'erasemode','xor')

    set(zh,'buttondownfcn','zpgui(''zeroclick'')',...
        'markersize',14,'linewidth',2)
    set(ph,'buttondownfcn','zpgui(''poleclick'')',...
        'markersize',14,'linewidth',2)
    
case 'addzero'
    if length(zh)>0
        zh(end+1) =  line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''zeroclick'')',...
        'markersize',14,'linewidth',4,'marker','o','linestyle','none');
    else
        zh = line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''zeroclick'')',...
        'markersize',14,'linewidth',4,'marker','o','linestyle','none');
        set(findobj('string','Remove Zero'),'enable','on')
    end
    zpgui('recompute')
    
case 'removezero'
    delete(zh(end))
    zh(end)=[];
    if length(zh)==0
        set(gco,'enable','off')
    end
    zpgui('recompute')
    
case 'addpole'
    if length(ph)>0
        ph(end+1) =  line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''poleclick'')',...
        'markersize',14,'linewidth',4,'marker','x');
    else
        ph =  line(.5,0,'parent',ax1,'buttondownfcn','zpgui(''poleclick'')',...
        'markersize',14,'linewidth',4,'marker','x');
        set(findobj('string','Remove Pole'),'enable','on')
    end
    zpgui('recompute')
    
case 'removepole'
    delete(ph(end))
    ph(end)=[];
    if length(ph)==0
        set(gco,'enable','off')
    end
    zpgui('recompute')

case 'zeroclick'
    
    set(gcf,'userdata','')
    set(gcf,'windowbuttonmotionfcn','set(gcf,''userdata'',''motion'')')
    set(gcf,'windowbuttonupfcn','set(gcf,''userdata'',''up'')')
    
    ind = find(zh==gco);
    set(zh(ind),'erasemode','xor')
    set(Lresp,'erasemode','xor')
    pair = length(get(zh(ind),'xdata'))==2;
    done = 0;
    
    while ~done
        waitfor(gcf,'userdata')
        switch get(gcf,'userdata')
        case 'motion'
            pt = get(ax1,'currentpoint');
            pt = pt(1,1:2);
            
            if pair
                set(zh(ind),'xdata',[pt(1) pt(1)],'ydata',[pt(2) -pt(2)])
            else
                set(zh(ind),'xdata',pt(1),'ydata',pt(2))
            end
            
            zpgui('recompute')
        case 'up'
            done = 1;
        end
        set(gcf,'userdata','')
    end
    set(gcf,'windowbuttonmotionfcn','')
    set(gcf,'windowbuttonupfcn','')
    set(zh(ind),'erasemode','normal')
    set(Lresp,'erasemode','normal')
    set(ax2,'ylimmode','auto')
    ylim = get(ax2,'ylim');
    Y = get(Lresp,'ydata');
    if ylim(1)>min(Y)-3 | ylim(2)<max(Y)+3
        set(ax2,'ylim',[min(Y)-3 max(Y)+3])
    end

case 'poleclick'
    
    set(gcf,'userdata','')
    set(gcf,'windowbuttonmotionfcn','set(gcf,''userdata'',''motion'')')
    set(gcf,'windowbuttonupfcn','set(gcf,''userdata'',''up'')')
    
    ind = find(ph==gco);
    set(ph(ind),'erasemode','xor')
    set(Lresp,'erasemode','xor')
    pair = length(get(ph(ind),'xdata'))==2;
    done = 0;
    
    while ~done
        waitfor(gcf,'userdata')
        switch get(gcf,'userdata')
        case 'motion'
            pt = get(ax1,'currentpoint');
            pt = pt(1,1:2);
            
            if pair
                set(ph(ind),'xdata',[pt(1) pt(1)],'ydata',[pt(2) -pt(2)])
            else
                set(ph(ind),'xdata',pt(1),'ydata',pt(2))
            end
            
            zpgui('recompute')
        case 'up'
            done = 1;
        end
        set(gcf,'userdata','')
    end
    set(gcf,'windowbuttonmotionfcn','')
    set(gcf,'windowbuttonupfcn','')
    set(ph(ind),'erasemode','normal')
    set(Lresp,'erasemode','normal')
    set(ax2,'ylimmode','auto')
    Y = get(Lresp,'ydata');
    ylim = get(ax2,'ylim');
    if ylim(1)>min(Y)-3 | ylim(2)<max(Y)+3
        set(ax2,'ylim',[min(Y)-3 max(Y)+3])
    end
    
case 'recompute'

    z = [];
    p = [];
    b = 1;
    a = 1;
    for i=1:length(zh)
        zx = get(zh(i),'xdata');
        zy = get(zh(i),'ydata');
        if length(zx)==1
            b = conv(b,[1 -(zx+sqrt(-1)*zy)]);
        else
            b = conv(b,[1 -2*zx(1) zx(1).^2+zy(1).^2]);
        end
        z = [z zx+sqrt(-1)*zy];
    end
    for i=1:length(ph)
        px = get(ph(i),'xdata');
        py = get(ph(i),'ydata');
        if length(px)==1
            a = conv(a,[1 -(px+sqrt(-1)*py)]);
        else
            a = conv(a,[1 -2*px(1) px(1).^2+py(1).^2]);
        end
        p = [p px+sqrt(-1)*py];
    end
    
    Y = fft(b,Nfft)./fft(a,Nfft);
    Y = Y/max(abs(Y));
    set(Lresp,'ydata',20*log10(fftshift(abs(Y))))

end

