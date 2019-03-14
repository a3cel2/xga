
/*  $Id: loader.js
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * File Description: Dynamic scripts loader and synchronizer
 */

timeStamp = new Date().getTime();
(function(){
    var isIE = (window.navigator.userAgent.indexOf('Trident/') >= 0);
    var refNode = document.scripts[document.scripts.length - 1];
    var srcjs = refNode.src.split('/');
    var viewer = ({sviewer: 'SeqView', treeviewer: 'TreeView', multialign: 'MultiAlignView'})[srcjs.pop().replace('.js', '')];
    if (window[viewer]) return;
    var dl = document.location;
    var externJS = (srcjs[2] != dl.hostname ? '-extern.js' : '.js');
    var dotJS = (srcjs[2].indexOf(/www|qa/) !== 0 || dl.href.indexOf('debug') > 0) ? '-debug.js' : '.js';
    
    var domNCBI = 'ncbi.nlm.nih.gov';
    var webNCBI = (dl.hostname.substr(-domNCBI.length) == domNCBI ? dl.protocol : 'https:') + '//'
        + (srcjs[2].indexOf(domNCBI) == -1 || srcjs[2].indexOf('blast.' + domNCBI) >= 0
        ? ('www.' + domNCBI) : srcjs[2]) + '/';

    var ExtJSver = document.location.href.indexOf('extjs621') > 0 ? '6.2.1' : '6.0.1';
    if (!document.querySelector('meta[name=viewport]')) {
        var meta = document.createElement('meta');
        meta.setAttribute('name', 'viewport');
        meta.setAttribute('content', 'width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no');
        document.getElementsByTagName('head')[0].appendChild(meta);
    }

    var appPath = ''
    for (var i = 3; i < srcjs.length - 1; i++)
        appPath += srcjs[i] + (srcjs[i] ? '/' : '');

    var extPath = webNCBI + 'core/extjs/ext-' + ExtJSver + '/build/';

    var thirdPartyScripts = {
        Ext: {tag: 'script', attr: {type: 'text/javascript', src: extPath + 'ext-all' + dotJS}},
        gapi: {tag: 'script', attr: {type: 'text/javascript', src: 'https://apis.google.com/js/client.js'}},
        jQuery: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + 'core/jquery/jquery-1.7.2.js'}},
        UUD: {tag: 'script',  attr: {type: 'text/javascript', src: webNCBI +'projects/genome/uud/js/uud' + externJS}},
        TMS: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + 'projects/genome/trackmgr/0.7/js/tms' + externJS}},
        ncbi: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + 'portal/portal3rc.fcgi/rlib/js/InstrumentOmnitureBaseJS/InstrumentNCBIBaseJS/InstrumentPageStarterJS.js'}},
        ExtCSS: {tag: 'link', attr: {rel: 'stylesheet', type: 'text/css', href: extPath + 'classic/theme-gray/resources/theme-gray-all.css'}},
        CSS: {tag: 'link', attr: {rel: 'stylesheet', type: 'text/css', href: webNCBI + appPath + 'css/style.css'}}};

    window[viewer + 'OnReady'] = function(callback, scope) {
        if (typeof window[viewer] === 'undefined')
            setTimeout(function() { window[viewer + 'OnReady'](callback, scope); }, 100);
        else {
            Ext.ariaWarn = Ext.emptyFn;
            Ext.onReady(callback, scope);
        }
    }

    var insertElement = function(resName, callback) {
        var res = thirdPartyScripts[resName];
        var deferCallback = !callback ? function(){} : function() {
            if (dotJS.length > 3) console.log(resName);
            if (!window[resName]) setTimeout(deferCallback, 30);
            else callback();
        }

        if (res.tag == 'link') {
            for (var p in document.styleSheets)
                if (res.attr.href == document.styleSheets[p].href) { res = false; break; }
        } else {
            if (typeof window[resName] !== 'undefined') res = false;
            else {
                for (var p in document.scripts)
                    if (res.attr.src == document.scripts[p].src) { res = false; break;}
            }
        }
        if (!res) {
           deferCallback();
           return;
        }

        var el = document.createElement(res.tag);
        for (var p in res.attr) el[p] = res.attr[p];
        if (typeof callback === 'function') {
            el.onload = el.onreadystatechange = function() {
                if (!this.readyState || this.readyState === "loaded" || this.readyState === "complete") {
                    if (isIE)
                        try {
                            var tmp = document.body.appendChild(document.createElement('div'));
                            tmp.innerText = ' ';
                            document.body.removeChild(tmp);
                        } catch(e){};
                    this.onreadystatechange = this.onload = null;
                    deferCallback();
                }
            }
        }
        refNode.parentNode.insertBefore(el, refNode);
    }
    insertElement('ExtCSS');
    insertElement('CSS');
    if (refNode.getAttribute('extjs') == 'skip') window.Ext = window.Ext || {};

    var counter = 100;

    insertElement('Ext', function finalLoad() {
        if ((typeof Ext === 'undefined' || !Ext.onReady) && counter--) {
            setTimeout(finalLoad, 100);
            return;
        }
        if (Ext.getVersion && Ext.getVersion().major >= 5) {
            Ext.enableAriaButtons = false;
            Ext.scopeCss = true;
            Ext.namespace(viewer);
            Ext.apply(window[viewer], {webNCBI: webNCBI, base_url: webNCBI + appPath, standalone: dl.pathname.indexOf(appPath) == 1});
            window['init' + viewer]();
        } else {
            alert('Current version of the Application works only in ExtJS ver. 6.+ environment!');
            return;
        }
     
        window[viewer].loadGAPI = function(url, callback) {
           insertElement('gapi', function(){ window[viewer].makeTinyURL(url, callback); });
        }
        if (Ext.isIE9) thirdPartyScripts['UUD'].attr.src = thirdPartyScripts['UUD'].attr.src.replace(externJS, '-extern.js');
        if (document.location.hostname.indexOf(domNCBI) == -1) {
            window[viewer].origHostname = document.location.hostname;
            __ncbi_stat_url = 'https://www.ncbi.nlm.nih.gov/stat';
            ncbi_signinUrl = 'https://www.ncbi.nlm.nih.gov/portal/signin.fcgi?cmd=Nop';
            ncbi_pingWithImage = true;    
        } else {
            if (webNCBI.indexOf('https://') != 0) {
               __ncbi_stat_url = 'https://dev.ncbi.nlm.nih.gov/stat';
               ncbi_pingWithImage = true;
            }
        }
        insertElement('jQuery', function(){
            insertElement('UUD', function(){ if (viewer == 'SeqView') insertElement('TMS'); });
        });
        if (!Ext.removeNodeOrig) {
            Ext.removeNodeOrig = Ext.removeNode;
            Ext.removeNode = function(n) { if (n) Ext.removeNodeOrig(n); }
        }

        window[viewer].extPath = extPath;
        window[viewer].extImagesPath = extPath + 'classic/theme-gray/resources/images/';
        window[viewer + 'OnReady'](function(){
            insertElement('ncbi');
            window[viewer].loadApp(refNode);
        });
    });
    
    if (window.NCBIGBUtils && NCBIGBUtils.makeTinyURL) return;
    
    NCBIGBUtils = window.NCBIGBUtils || {};
    NCBIGBUtils.makeTinyURL = function(url, callback) {
        var urlService = function() {
            var request = gapi.client.urlshortener.url.insert({'resource': {'longUrl': url}});
            request.execute(function(res) {
                res.id = (res.id || '').replace(/http:/, 'https:');
                callback(res);
            });
        }
        if (typeof gapi === 'undefined') { // load gapi
            insertElement('gapi', function(){ NCBIGBUtils.makeTinyURL(url, callback); });
            return; 
        }
        if (!gapi.client) {
            if (!gapi.load) callback({id: false}); //not loaded
            else // wait while loading
                setTimeout(function(){NCBIGBUtils.makeTinyURL(url, callback);}, 500);
            return;
        }
        if (!gapi.client.urlshortener) {
            var googleAPIkey = 'AIzaSyC998TlAgTLRQd91sRa2NJDDE8ieXrr5XQ';
            gapi.client.setApiKey(googleAPIkey);
            gapi.client.load('urlshortener', 'v1', urlService);
        } else urlService();
    }
})();




/*  $Id: globals.js 37949 2017-03-07 21:03:19Z borodine $
 * Authors:  Maxim Didenko
 * File Description:
 *
 */

function initSeqView()
{

SeqView.doInitPing = true;
SeqView.loadApp = function (refNode) {
    var items = Ext.query('div[class=SeqViewerApp]');
    for (var i = 0; i < items.length; i++) {
        var item = items[i];
        var id = item.id || Ext.id();
        if (!item.hasAttribute('data-autoload') && refNode.id != 'autoload') continue; 
        if (SeqView.App.findAppByDivId(id)) return; 
        item.id = id;
        var app = new SeqView.App(id);
        app.load();
    }
}

if (Ext.isIE) {
    var db = document.querySelector('meta[name=ncbi_db]');
    SeqView.geneIE = (db && db.content && (db.content == 'clone' || db.content == 'gene'));
}

Ext.define('SeqView.SearchField', {
    extend: 'Ext.form.field.Text',
    width:180,
    triggers: { 
        clear: {
            cls: 'x-form-clear-trigger',
            hidden: true,
            handler: function(tb, trg) { this.clickTrigger(trg.id); }
        },
        search: {
            cls: 'x-form-search-trigger',
            handler: function(tb, trg) { this.clickTrigger(trg.id); }
        }
    },
    initComponent: function(){
        this.callParent(arguments);
        this.on('specialkey', function(f, e){
            if(e.getKey() == e.ENTER){
                this.clickTrigger(null, {id: 'search'});
            }
        }, this);
    },
    

    clickTrigger: function(trgid) {
        if (trgid == 'clear') {
            this.triggers.clear.hide();
            this.setValue('');
        } else this.triggers.clear.show();
        var bp = this.store.baseParams;
        bp = bp || {};
        bp.query = this.getRawValue();
        bp.regexp = bp.query ? globStringToRegex(bp.query) : '';
        if (bp.filterBy) {
            this.store.filterBy(bp.filterBy, this.store);
        } else this.store.reload({params:{start: 0}});
        
    }
});

if (navigator.maxTouchPoints > 1 || 'ontouchstart' in document) {
    SeqView.useTouch = true;
}

SeqView.Cookies = function(config) {
    this.path = "/";
    this.expires = new Date(new Date().getTime()+(1000*60*60*24*7)); //7 days
    this.domain = ".nih.gov"; // this is essential for proper cookie functioning
    this.secure = false;
    Ext.apply(this, config);

    var readCookies = function() {
        var cookies = {};
        var c = document.cookie + ";";
        var re = /\s?(.*?)=(.*?);/g;
        var matches;
        while((matches = re.exec(c)) != null){
            var name = matches[1];
            // fix incorect cookie name for IE8 (like 'WebCubbyUser; sv-user-data')
            var a = name.split(';');
            name = a[a.length-1].trim();
            var value = matches[2];
            cookies[name] = value;
        }
        return cookies;
    }
    this.state = readCookies();
};

SeqView.Cookies.MaxDate = new Date(new Date('10000/01/01 GMT').getTime()-100);
SeqView.Cookies.UserTracksCookieNameBase = 'sv-usertracks-key';
SeqView.Cookies.UserDataCookieNameBase = 'sv-userdata-key';
SeqView.Cookies.AppDataCookieNameBase = 'sv-data-key';
SeqView.Cookies.UserTracksCookieName = SeqView.Cookies.UserTracksCookieNameBase;
SeqView.Cookies.UserDataCookieName = SeqView.Cookies.UserDataCookieNameBase;
SeqView.Cookies.AppDataCookieName = SeqView.Cookies.AppDataCookieNameBase;

SeqView.Cookies.prototype = {
    get : function(name, defaultValue){
        return typeof this.state[name] == "undefined" ?
            defaultValue : this.state[name];
    },

    set : function(name, value){
        if(typeof value == "undefined" || value === null){
            this.clear(name);
            return;
        }
        this.setCookie(name, value);
        this.state[name] = value;
    },

    clear : function(name){
        delete this.state[name];
        this.clearCookie(name);
    },

    // private
    setCookie : function(name, value){
        document.cookie = name + "=" + value +
           ((this.expires == null) ? "" : ("; expires=" + this.expires.toGMTString())) +
           ((this.path == null) ? "" : ("; path=" + this.path)) +
           ((this.domain == null) ? "" : ("; domain=" + this.domain)) +
           ((this.secure == true) ? "; secure" : "");
    },

    // private
    clearCookie : function(name){
        document.cookie = name + "=del; expires=Thu, 01-Jan-1970 00:00:01 GMT" +
           ((this.path == null) ? "" : ("; path=" + this.path)) +
           ((this.domain == null) ? "" : ("; domain=" + this.domain)) +
           ((this.secure == true) ? "; secure" : "");
    }
};
SeqView.SessionData = new SeqView.Cookies({expires:null, secure:Ext.isSecure});
SeqView.UserTracks = new SeqView.Cookies({expires:SeqView.Cookies.MaxDate, secure:Ext.isSecure});
SeqView.trackSets = {tmsSets: [], usrSets: [], trackLists: []};

SeqView.configureTrackSets = function(app) {
    TMS.TrackSets.TrackSetService.SetAssembly(app.m_AssmContext);
    TMS.TrackSets.GetDefaultTrackSet('SViewer', app.m_AppContext, app.m_AssmContext)
        .done(function() {
            app.m_defaultTrackSet = this.GetTracks();
        })
        .fail(function() {
            console.log('Failed to get default track set');
        });
    SeqView.requestTrackSets = function(callback) {
        TMS.TrackSets.TrackSetService.GetTracksets(true)
           .done(callback)
           .fail(function() { callback([]);}
        );
    }
};

SeqView.AreaFlags = {
    Link:           (1 << 0), ///< a direct link stored in m_Action
    CheckBox:       (1 << 1), ///< a toggle button/icon
    NoSelection:    (1 << 2), ///< the object can't be selected
    ClientSelection:(1 << 3), ///< the selection can be done on client
    NoHighlight:    (1 << 4), ///< on highlighting on mouse over
    NoTooltip:      (1 << 5), ///< do not request and show tooltip
    TooltipEmbedded:(1 << 6), ///< tooltip embedded
    Track:          (1 << 7), ///< track title bar
    Ruler:          (1 << 8),  ///< ruler bar
    Editable:       (1 << 9),  ///< editable area
    NoPin:          (1 << 10), ///< not pinnable
    IgnoreConflict: (1 << 11), ///< feature can be ignored (isca browser feature editing only)
    Sequence:       (1 << 12), ///< sequence track flag
    Comment:        (1 << 13), ///< render a label/comment on client side
    DrawBackground: (1 << 14), ///< highlight background for this area
    Dirty:          (1 << 16), ///< dirty flag
    NoNavigation:   (1 << 17), ///< no havigation buttons on title bar
    Legend:         (1 << 18)  ///< legends for graph_overlay tracks
};


SeqView.getHelpURL = function() {
    return 'https://www.ncbi.nlm.nih.gov/tools/sviewer/';  
};

SeqView.showHelpDlg = function(extra_params) {
    window.open(SeqView.getHelpURL() + (extra_params ? extra_params : ''));
};

SeqView.showAboutMessage = function() {
    var msg = '<p>NCBI Graphical Sequence Viewer - graphical display for the Nucleotide and Protein sequences.</p>Version '
    + SeqView.getVersionString() + '<br><br><span style="font-size:10px;">' + SeqView.version.rev + '</span><br>';
    msg += '<br>ExtJS version: ' + Ext.getVersion('core').version + '<br>';
    msg += '<br><a href=\"' + SeqView.base_url + 'info.html' + '\" target=\"_blank\" style=\"color:blue\">CGI binaries Info</a>';
    Ext.Msg.show({
        title: 'Sequence Viewer',
        msg: msg,
        maxWidth: '320',
        minWidth: '300',
        buttons: Ext.Msg.OK
    });
}

SeqView.getVersion = function() { return parseFloat(SeqView.version.num); }
SeqView.getVersionString = function() { return SeqView.version.str }

SeqView.decode = function(data) { return (typeof data === 'object') ? data : Ext.decode(data); }

SeqView.m_Apps = [];


SeqView.version = {rev: 'Revision: 38416, install date: 2017-05-08 11:10', str: '3.21.0', num:'3210'};
/* $Id: utils.js
 * File Description: Univirsal tools and objects
 */
if (!window.NCBIGBUtils || !window.NCBIGBUtils.ClearBrowserSelection) {
var utils = window.NCBIGBUtils = window.NCBIGBUtils || {};
// Decorate long numbers with commas
Number.prototype.commify = function() {
    nStr = this + '';
    x = nStr.split('.');
    x1 = x[0];
    x2 = x.length > 1 ? '.' + x[1] : '';
    var rgx = /(\d+)(\d{3})/;
    while (rgx.test(x1)) {
        x1 = x1.replace(rgx, '$1' + ',' + '$2');
    }
    return x1 + x2;
};

Number.prototype.shorten = function() {
    var value = this.valueOf();
    var negv = ( value < 0 );
    if( negv ) value = -value;
    var Suffixes = [ '', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y' ];
    var sfx = '';
    for( var i = 0; i < Suffixes.length; i++ ){
        if( value < 1000 ){
            sfx = Suffixes[i];
            break;
        }
        value /= 1000;
    }
    return ( negv ? '-' : '' ) + Number( value ).toFixed( value < 10 && sfx != '' ? 1 : 0 ) + sfx;
};

/* Cross-Browser Split 1.0.1
(c) Steven Levithan <stevenlevithan.com>; MIT License
An ECMA-compliant, uniform cross-browser split method */

var cbSplit;

// avoid running twice, which would break `cbSplit._nativeSplit`'s reference to the native `split`
if (!cbSplit) {

cbSplit = function (str, separator, limit) {
    // if `separator` is not a regex, use the native `split`
    if (Object.prototype.toString.call(separator) !== "[object RegExp]") {
        return cbSplit._nativeSplit.call(str, separator, limit);
    }

    var output = [],
        lastLastIndex = 0,
        flags = (separator.ignoreCase ? "i" : "") +
                (separator.multiline  ? "m" : "") +
                (separator.sticky     ? "y" : ""),
        separator = RegExp(separator.source, flags + "g"), // make `global` and avoid `lastIndex` issues by working with a copy
        separator2, match, lastIndex, lastLength;

    str = str + ""; // type conversion
    if (!cbSplit._compliantExecNpcg) {
        separator2 = RegExp("^" + separator.source + "$(?!\\s)", flags); // doesn't need /g or /y, but they don't hurt
    }

    /* behavior for `limit`: if it's...
    - `undefined`: no limit.
    - `NaN` or zero: return an empty array.
    - a positive number: use `Math.floor(limit)`.
    - a negative number: no limit.
    - other: type-convert, then use the above rules. */
    if (limit === undefined || +limit < 0) {
        limit = Infinity;
    } else {
        limit = Math.floor(+limit);
        if (!limit) {
            return [];
        }
    }

    while (match = separator.exec(str)) {
        lastIndex = match.index + match[0].length; // `separator.lastIndex` is not reliable cross-browser

        if (lastIndex > lastLastIndex) {
            output.push(str.slice(lastLastIndex, match.index));

            // fix browsers whose `exec` methods don't consistently return `undefined` for nonparticipating capturing groups
            if (!cbSplit._compliantExecNpcg && match.length > 1) {
                match[0].replace(separator2, function () {
                    for (var i = 1; i < arguments.length - 2; i++) {
                        if (arguments[i] === undefined) {
                            match[i] = undefined;
                        }
                    }
                });
            }

            if (match.length > 1 && match.index < str.length) {
                Array.prototype.push.apply(output, match.slice(1));
            }

            lastLength = match[0].length;
            lastLastIndex = lastIndex;

            if (output.length >= limit) {
                break;
            }
        }

        if (separator.lastIndex === match.index) {
            separator.lastIndex++; // avoid an infinite loop
        }
    }

    if (lastLastIndex === str.length) {
        if (lastLength || !separator.test("")) {
            output.push("");
        }
    } else {
        output.push(str.slice(lastLastIndex));
    }

    return output.length > limit ? output.slice(0, limit) : output;
};
cbSplit._compliantExecNpcg = /()??/.exec("")[1] === undefined; // NPCG: nonparticipating capturing group
cbSplit._nativeSplit = String.prototype.split;
} // end `if (!cbSplit)`

// for convenience... interferes with Ext.js interpretation of split, sorry
//String.prototype.split = function (separator, limit) {
//    return cbSplit(this, separator, limit);
//};
// End of Cross-Browser Split

String.prototype.trimToPix = function(length) {
    var tmp = trimmed = utils.sanitize(this);
    if (tmp.visualLength() > length)  {
        trimmed += "...";
        while (trimmed.visualLength() > length)  {
            tmp = tmp.substring(0, tmp.length-1);
            trimmed = tmp + "...";
        }
    }
    return trimmed;
};
String.prototype.trim = function () {
    return this.replace(/^\s*/, "").replace(/\s*$/, "");
};

utils.ClearBrowserSelection = function() {
    var sel;
    try {
        if (document.selection && document.selection.empty) document.selection.empty();
        else 
            if (window.getSelection) {
                sel = window.getSelection();
                if (sel && sel.removeAllRanges) sel.collapseToEnd();
            }
    } catch(e) {}
};

utils.isNumeric = function(str) { return /^-?[0-9]+(\.[0-9]*)?[km]?$/i.test(str); };

utils.stringToNum = function(pos_str) {
    if (!pos_str) return;
    pos_str = pos_str.replace(/[, ]/g,'');
    if (pos_str.length < 1) return;
    var multiplier = 1;
    var last_char = pos_str.charAt(pos_str.length - 1).toUpperCase();
    if (last_char == 'K' || last_char == 'M') {
        pos_str = pos_str.substr(0, pos_str.length - 1);
        if (last_char == 'K') {
           multiplier = 1000;
        } else {
           multiplier = 1000000;
        }
    }
    var dec_part = 0;
    if (multiplier > 1) {
        var dec_pos = pos_str.indexOf('.');
        if (dec_pos != -1) {
            dec_part = Math.floor(parseFloat(pos_str.substr(dec_pos)) * multiplier);
            pos_str = pos_str.substr(0, dec_pos);
        }
    }
    return multiplier * parseInt(pos_str) + dec_part;
};

// Escape handles symbols innocuous from the point of view of HTML/HTTP
// but used internally in names, which are ususally provided by the user
// and as such are out of control. Encoding schema is according to SV-591
// and SV-1379:
// Prepend \ | , : [ ] with backslash, encode & ; # % as \hex-code, encode
// ' " = and space using standard %hex notation.
// Added a parameter 'symbs' allowing to process only a specific set of symbols
utils.escapeName = function(s, symbs) {
    var res = "";
    var parts = [];
    if (!symbs) parts = cbSplit(s, /([\]\[\\\|\'\"= ,:;&#%])/);
    else parts = cbSplit(s, symbs);
    for (var len = parts.length, i = 0; i < len;) {
        res += parts[i++];
        var sym = parts[i++];
        if (!sym) break;
        if (/[\]\[\\\|,:]/.test(sym))
            res += "\\" + sym;
        else {
            if (/['"= ]/.test(sym))
                res += "%";
            else
                res += "\\"; // ; & # %
            res += sym.charCodeAt(0).toString(16);
        }
    }
    return res;
};

utils.escapeTrackName = function(s, symbs) {
    var res = "";
    if (!symbs) symbs = /([\]\[\\\|,:=&;"#%])/;
    var parts = cbSplit(s, symbs);
    for (var len = parts.length, i = 0; i < len;) {
        res += parts[i++];
        var sym = parts[i++];
        if (!sym) break;
        res += "\\";
        if (/[=&\\;"#%]/.test(sym)) {
            res += sym.charCodeAt(0).toString(16);
        } else {
            res += sym;
        }
    }
    return res;
};

utils.unescapeName = function(s) {
    var parts = unescape(s).split("\\");
    var res = parts[0];
    for (var len = parts.length, i = 1; i < len; i++) {
        var part = parts[i];
        if (part.length == 0)
            res += "\\" + parts[++i];
        else if (/[\]\[\\\|,:]/.test(part.charAt(0)))
            res += part;
        else if (/^[0-9a-f]{2}.*/.test(part)) {
            res += String.fromCharCode(parseInt(part.slice(0,2), 16)) + part.slice(2);
        }
    }
    return res;
};


utils.sanitize = function(s) {
    return s.replace(/&/g, "&amp;").replace(/</g, "&lt;")
            .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
};

String.prototype.visualLength = function() {
   var ruler = document.getElementById('string_ruler_unit');
   ruler.innerHTML = this;
   var ret = ruler.offsetWidth;
   ruler.innerHTML = '';
   return ret;
};
} //if (window.NCBIGBUtils)
/*  $Id: tooltip.js
 * File Description: Ext.ToolTip extention
*/
if (!window['NCBIGBObject.ToolTip']) {
Ext.tip.QuickTipManager.init(true);// , {showDelay: 1000, mouseOffset: [-15, 15]});

if (Ext.supports.Touch && Ext.versions.core.version.indexOf('6.0.') == 0) {
    Ext.define('Ext.tip.ToolTip', {
        override: 'Ext.tip.ToolTip',
        setTarget: function(target) {
            var me = this,
                t = Ext.get(target),
                tg;
            if (me.target) {
                tg = Ext.get(me.target);
                if (Ext.supports.Touch) me.mun(tg, 'tap', me.onTargetOver, me);
                me.mun(tg, {
                    mouseover: me.onTargetOver,
                    mouseout: me.onTargetOut,
                    mousemove: me.onMouseMove,
                    scope: me
                });
            }
            me.target = t;
            if (t) {
                if (Ext.supports.Touch)  me.mon(t, { tap: me.onTargetOver, scope: me });
                me.mon(t, {
                    mouseover: me.onTargetOver,
                    mouseout: me.onTargetOut,
                    mousemove: me.onMouseMove,
                    scope: me
                });
            }
            if (me.anchor) {
                me.anchorTarget = me.target;
            }
        }
    });
}

Ext.define('NCBIGBObject.ToolTip', {
    extend: 'Ext.tip.ToolTip',
    hideDelay: 300,
    showDelay: 1000,
    draggable: true,
    pinned: false,
    autoScroll: false,
    adjustWidth: false,
    pinCallbackFn: Ext.emptyFn,
    header: {
        baseCls: 'x-panel-header',
        padding: 0,
        titlePosition: 99
    },
    listeners: {
        beforehide: function() {
            return !(this.insideTT || this.pinned);
        }
    },
    initComponent: function() {
        if (this.pinnable) delete this.pinnable;
        this.callParent(arguments);
        if (this.isPinnable())
            this.addTool([{type: 'unpin', scope: this, hidden: false, callback: this.toggelPin}]);
     },
/*
     afterRender: function() {
         this.callParent(arguments);
//         this.body.setStyle('position', 'absolute');
    },*/
    isPinnable: function() { return this.pinnable != false; },

    doWidthAdjusting: function() {
        if (!this.adjustWidth || !this.isVisible()) return; 
        var maxW = 490;
            minW = this.getWidth();
        if (maxW < minW || maxW - minW < 50 || minW <= this.minWidth) return;

        var delta = 0,
            curH = this.getHeight(),
            minH = this.setWidth(maxW).getHeight();

        while (minH < curH) {
            delta = maxW - minW;
            delta = delta >> 1;
            curH = this.setWidth(minW + delta).getHeight();
            if (minH < curH++) minW += delta;
            else maxW -= delta;

            if (delta < 3) break;
        }
        this.setWidth(minW + delta);
    },

    onMouseOver: function(e) {
        e.stopEvent();
        this.insideTT = true;
    },

    onMouseOut: function(e) {
        e.stopEvent();
        this.insideTT = e.within(this.el, true, true);
        if (!(this.insideTT || this.pinned)) this.hide();
    },

    onShow: function() {
        this.callParent(arguments);
        this.getEl().on('mouseover', this.onMouseOver, this);
        this.getEl().on('mouseout', this.onMouseOut, this);
        this.doWidthAdjusting();
    },

    toggelPin: function(o, tool) {
        if (!this.pinned) {
            this.pinCallbackFn();
            tool.setType('pin');
            if (!this.isVisible()) this.show();
        } else {
            tool.setType('unpin');
            NCBIGBObject.ToolTip.superclass.onHide.call(this);
        }
        this.autoHide = !(this.pinned = !this.pinned);
    },
    update: function(html) {
        if (this.adjustWidth) delete this.width; // = 'auto';
        NCBIGBObject.ToolTip.superclass.update.call(this, html);
        this.updateLayout();
        this.doWidthAdjusting();
    }
});
}

/*  $Id: locator.js 38012 2017-03-14 20:19:18Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko
 * File Description: Locator for the overview panel & alignment panel
 */

/* locator */
SeqView.Locator = function(view, color, resizable) {
    this.m_View = view; // view 
    this.m_panorama = Ext.get(view.m_App.m_Views[0].m_DivId);

    var tmpl_string_pan = '<div id="pan_scroller_{idx}" class="pan-bar" style="background-color:#{color};position:absolute;top:0px;left:0px;"></div>';
    var tmpl_string = '<div id="seq_g{idx}" style="position:absolute;top:17px;left:0px;width:1px;">';
    tmpl_string += '<div class="locator_rect"></div>';
    if (resizable) {
        tmpl_string += '<div id="left_resizer_{idx}" class="left-resizer"></div><div id="right_resizer_{idx}" class="right-resizer"></div>';
    }
    tmpl_string += '</div>';
      
    var idx = this.m_View.m_Idx;     
    var tpl_pan = new Ext.Template(tmpl_string_pan);
    var tpl = new Ext.Template(tmpl_string);
    var panorama_height = this.getPanoramaHeight();
    this.element_pan = tpl_pan.insertFirst(this.m_panorama, {idx:idx, qtip:'View #', color:color, adjheight: panorama_height-19}, true);
    this.element = tpl.insertFirst(this.m_panorama, {idx:idx, qtip:'View #', color:color, adjheight: panorama_height-19}, true);
    this.element.setHeight(panorama_height - 2);

    this.element_pan.on({
        'mousedown'  : this.onMouseDown,
        'touchstart' : this.onMouseDown,
        'contextmenu':  this.onContextMenu,  
        scope: this
    });
    
    if (resizable) {  
        Ext.get("left_resizer_" + idx).on({
            'mousedown' : this.onMouseDown,
            'touchstart': this.onMouseDown,
            scope: this
        });
        Ext.get("right_resizer_" + idx).on({
            'mousedown' : this.onMouseDown,
            'touchstart': this.onMouseDown,
            scope: this
        });
    }
};

SeqView.Locator.prototype = {
    getElement:  function() { return this.element; },    
    setColor:  function(color) {
        Ext.get('pan_scroller_'+this.m_View.m_Idx).setStyle('background-color', '#' + color);
        Ext.get(this.m_View.m_DivId).parent().setStyle('borderLeft', '#' + color + ' 1px solid');
    },
    setHeight:  function(h) { 
        this.element.setHeight(h - 17);  
    },
    getLeft:  function(local) { return this.element.getLeft(local); },
    setLeft:  function(pos) { 
        var panorama_width = this.getPanoramaWidth();
        var bar_pos = pos;
        if(pos + 24 > panorama_width)
            bar_pos = panorama_width - 24;

        this.element_pan.setLeft(bar_pos);
        this.element.setLeft(pos); 
    },
    getWidth:  function(contentWidth) { return this.element.getWidth(contentWidth); },
    setWidth: function(width) { 
        this.element.setWidth(width); 
    }, 
    getRight:  function(local) { return this.element.getRight(local); },
    remove:  function() { 
        this.element_pan.remove(); 
        this.element.remove(); 
    },

    getPanoramaWidth: function() {
        return this.m_View.m_App.getPanoramaWidth();
    },
    getPanoramaHeight: function() {
        return this.m_View.m_App.getPanoramaHeight();
    },


//////////////////////////////////////////////////////////////////////////
// onMouseDown:

    onMouseDown: function(e) {
        if (this.m_ContextMenu) this.m_ContextMenu.destroy();
        if ((e.type == 'mousedown' && e.button) || this.m_View.m_App.m_Panorama.m_Loading) return;
        var elem = e.getTarget().id;
        this.m_XY = e.getXY();
        this.m_Action = 'Resize';         
        if (elem.indexOf('left_resizer') == 0) {
            this.m_moveHandler = function(delta) {
                var old_left = this.getLeft(true);
                var new_left = Math.max(0, old_left - delta);
                new_left = Math.min(this.getRight(true) - 2, new_left);
                var new_width = this.getWidth() + (old_left - new_left);
                if (this.m_View.checkLocatorWidth(new_width) ) {
                    this.setLeft(new_left);
                    this.setWidth(new_width);
                    this.m_Action = 'Resize';
                }
            }
        } else if (elem.indexOf('right_resizer') == 0) {
            this.m_moveHandler = function(delta) {
                var new_width = Math.max(2, this.getWidth() - delta);
                new_width = Math.min(this.getPanoramaWidth() - 2 - this.getLeft(true), new_width);
                if (this.m_View.checkLocatorWidth(new_width) ) {
                    this.setWidth(new_width);
                    this.m_Action = 'Resize';
                }
            }
        } else if (elem.indexOf('pan_scroller') == 0) {
            if (e.type != 'mousedown')
                this.m_deferredContext = Ext.defer(SeqView.fireEvent, 2000, this, ['contextmenu', elem]);
            this.m_moveHandler = function(delta) {
                var new_left = Math.max(0, this.getLeft(true) - delta);
                new_left = Math.min(this.getPanoramaWidth() - 2 - this.getWidth(), new_left);
                this.setLeft(new_left);
                this.m_Action = 'Drag';     
            }
        } else { return; }
            
        var locator = this;
        var onMove = function(e) {
            locator.onMouseMove(new Ext.event.Event(e ? e : window.event));
        }
        var onEnd = function(e) {
            locator.onMouseUp(new Ext.event.Event(e ? e : window.event));
        }
        if (e.button == 0) {
            this.m_DocMouseMove = document.onmousemove;
            this.m_DocMouseUp = document.onmouseup;
            document.onmousemove = onMove;
            document.onmouseup = onEnd;
        } else {
            this.m_DocTouchMove = document.ontouchmove;
            this.m_DocTouchUp = document.ontouchend;
            document.ontouchmove = onMove;
            document.ontouchend = onEnd;
        }                
        if (Ext.isIE) e.stopPropagation(); else e.stopEvent();
    },

    onMouseUp: function(e) {
        if (this.m_Action) {
            this.m_View.syncToLocator();
            this.m_View.pingClick('1-0-' + this.m_Action.charAt(0));
        }
        if (this.m_deferredContext) {   
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        this.m_moveHandler = function(){};
        this.m_Action = '';  
        this.m_XY = null;
        if (e.button == 0) {
            document.onmousemove = this.m_DocMouseMove;
            document.onmouseup = this.m_DocMouseUp;
        }
        else {
            document.ontouchmove = this.m_DocTouchMove;
            document.ontouchend = this.m_DocTouchUp;            
        }
        this.m_DocMouseMove = this.m_DocMouseUp = this.m_DocTouchMove = this.m_DocTouchUp = null;
        e.stopPropagation();
    },
    
    onMouseMove: function(e) {
        if (!this.m_XY) return;
        NCBIGBUtils.ClearBrowserSelection();
        var majicSensivity = 5; 
        var xy = e.getXY();
        if (e.browserEvent.changedTouches && this.m_deferredContext) {
            if (Math.abs(xy[0] - this.m_XY[0]) > majicSensivity
             || Math.abs(xy[1] - this.m_XY[1]) > majicSensivity)
            {   
                clearTimeout(this.m_deferredContext);
                this.m_deferredContext = 0;
            }
        }
        this.m_moveHandler(this.m_XY[0] - xy[0]);
        this.m_XY = xy; // save new values        
        e.stopPropagation();
    },
/*
    onMouseOut: function(e) {
        //console.log(e);
    },*/

    onContextMenu: function(e) {
        e.stopEvent();
        var menu = new Ext.menu.Menu();
        var elem = e.getTarget().id;
        if (!this.m_deferredContext) {   
            this.m_XY = e.getXY();
        } else {
            this.m_deferredContext = 0;
            this.m_ContextMenu = menu;
        } 
        var pan_scroller = Ext.get(elem);
        var idx = elem.split('_')[2];
        var seq_id = 'seq_g' + idx;
        menu.add([{ 
                text: 'Bring to front', handler:function() {             
                    var pan_holder = Ext.get('pan-holder');
                    var seq_g = Ext.get(seq_id);
                    seq_g.insertBefore(pan_holder);
                    pan_scroller.insertAfter(seq_g);
                    }}, {
                text: 'Sent to back', handler:function() {  
                    var seq_id = 'seq_g' + idx;
                    var seq_g = Ext.get(seq_id);
                    this.m_panorama.insertFirst(seq_g)
                    pan_scroller.insertAfter(seq_g);
                    }, scope: this }, {
                text: 'View color change',
                    menu:new Ext.menu.ColorPicker({listeners: {'select': function(cm, color) {
                        this.m_View.m_Color = color;
                        this.setColor(color);
                        this.m_View.m_App.reCreateReflections();
                    }, scope: this }})
                }
            ]);
        menu.showAt(this.m_XY);
    }
};


/*  $Id: reflection.js 35303 2016-04-21 16:39:21Z rudnev $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vlad Lebedev, Maxim Didenko
 *
 * File Description: Visible range reflections in the viewer
 *
 */
 
 
/* view reflections */

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Reflection 
/********************************************************************/

 SeqView.Reflection = function(this_view, show_view) {
    this.m_PrevXY = [];

    this.m_View = this_view; // view index
    this.m_ShowView = show_view; // view index
    
    var color = this.m_ShowView.m_Color;
    
    this.parent_id = this.m_View.m_DivId;
    this.parent_elem = Ext.get( this.parent_id );
    
    var vis_range = this.m_ShowView.toSeq();
    var qtip = ( vis_range[0] +1 ) + ' : ' + ( vis_range[1] +1 );
    
    var tpl = new Ext.Template(
	    '<div id="reflection_id_{idx}_{sidx}" data-qtip="{qtip}" class="reflection" style="z-index:10;background-color:#{color};"/>'
	);
    this.element = tpl.append(
		this.parent_elem, 
		{idx:this.m_View.m_Idx,sidx:this.m_ShowView.m_Idx, qtip:qtip, color:color}, 
		true
	);
    
	this.element.on({
        'mousedown' : this.onMouseDown,
        //'mouseup' :   this.onMouseUp,
        //'mousemove' : this.onMouseMove,
        scope:this
    });

 	this.update();
}


SeqView.Reflection.prototype = {

	update: function() {
	
		if( !this.m_View || !this.m_ShowView ) return;
	
		var show_range = this.m_ShowView.toSeq();
		// conversion in pixels
		var f = this.m_View.seq2Pix( show_range[0] ); 
		var t = this.m_View.seq2Pix( show_range[1] );

		if( this.m_View.getFlip() ){
			this.element.setLeft( Math.round( t +this.m_View.m_ScrollPix ) ); // mind the scroll offset!
			this.element.setWidth( Math.round( f -t ) );
		} else {
			this.element.setLeft( Math.round( f +this.m_View.m_ScrollPix ) ); // mind the scroll offset!
			this.element.setWidth( Math.round( t -f ) );
		}
		
		var tooltip_id = "reflection_id_" + this.m_View.m_Idx + "_" + this.m_ShowView.m_Idx;
		var tooltip = Ext.get( tooltip_id );

		if( tooltip ){
			var qtip = (show_range[0]+1) + '-' + (show_range[1]+1);
			var qtsetup = {};
			qtsetup["data-qtip"] = qtip;
			tooltip.set( qtsetup );
		}	
	},

    movePix: function(delta) {
        var new_left = this.element.getLeft(true) - delta;
        this.element.setLeft(new_left);
    },
    
    remove: function() { this.element.remove(); },
    
    scrollPix: function(view, delta) {
        if (this.m_View.m_Idx == view.m_Idx) { // scroll as one (locked)
            var new_left = this.element.getLeft(true) - delta;
            this.element.setLeft(new_left);
        } else if (this.m_ShowView.m_Idx == view.m_Idx) {
            var vis_range = view.toSeq();
            var f = this.m_View.seq2Pix(vis_range[0]); // in pixels
            this.element.setLeft( f + this.m_View.m_ScrollPix ); // mind the scroll offset!
        }
    },
    
//////////////////////////////////////////////////////////////////////////
// onMouseDown:

    onMouseDown: function(e) {
        if (e.button === 0) {
            this.m_PrevXY = e.getXY(); 
            this.element.setStyle('cursor', 'move');
            this.m_DocMouseMove = document.onmousemove;
            this.m_DocMouseUp = document.onmouseup;
            var refl = this;
            document.onmousemove = function(e) { refl.onMouseMove(new Ext.event.Event(e ? e : window.event)); };
            document.onmouseup = function(e) { refl.onMouseUp(new Ext.event.Event(e ? e : window.event)); };
        }
    },

//////////////////////////////////////////////////////////////////////////
// onMouseUp:

    onMouseUp: function(e) {
        if (!this.m_PrevXY) { return; }
        this.element.setStyle('cursor', 'pointer');
        var flip = this.m_View.getFlip();
        
        var pix_from = Math.abs(this.m_View.m_ScrollPix) + this.element.getLeft(true);
        var pix_to   = pix_from + this.element.getWidth();
              
        var from_seq = this.m_View.pix2Seq(flip ? pix_to : pix_from);
        var to_seq   = this.m_View.pix2Seq(flip ? pix_from : pix_to);
        
        this.m_ShowView.startImageLoading(from_seq, to_seq-from_seq+1);
        
        document.onmousemove = this.m_DocMouseMove;
        document.onmouseup = this.m_DocMouseUp;
        this.m_DocMouseMove = null;
        this.m_DocMouseUp = null;
        this.m_PrevXY = null;
    },

//////////////////////////////////////////////////////////////////////////
// onMouseMove:

    onMouseMove: function(e) {
        if (!this.m_PrevXY) { return; }
        var xy = e.getXY();
        var delta_x = this.m_PrevXY[0] - xy[0];
        var new_left = this.element.getLeft(true) - delta_x;
        this.element.setLeft(new_left);
        this.m_PrevXY = xy;
    }
    
};

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.ReflectionCont 
/********************************************************************/

SeqView.ReflectionCont = function(app) {
    this.m_App = app;
    this.m_AllReflections = [];
    
};

SeqView.ReflectionCont.prototype = {
    
    reCreate: function() {
        this.deleteAll();

        this.m_App.forEachView(function(this_view) {
            if (this_view.isGraphic()) {
                this.m_App.forEachView(function(show_view) { // create new ones
                    if (show_view.isGraphic() &&  show_view.m_Idx != this_view.m_Idx) {
                        this.m_AllReflections.push(new SeqView.Reflection(this_view, show_view));
                    }
                }, this);
                if (this.m_App.m_TextView) { // sequence text  view shown
                    this.m_AllReflections.push(new SeqView.Reflection(this_view, this.m_App.m_TextView));
                }
            }
        },this);
    },

	updateAll: function() {
        for( var i = 0; i < this.m_AllReflections.length; i++ ){ 
            this.m_AllReflections[i].update();
        }
	},
   
    deleteAll: function() {
        for (var i = 0; i != this.m_AllReflections.length; i++) { 
            this.m_AllReflections[i].remove();
        }
        this.m_AllReflections = [];
    },
    
    scrollPix: function(view, delta) {
        for (var i = 0; i != this.m_AllReflections.length; i++) {
            var reflection = this.m_AllReflections[i];
            reflection.scrollPix(view, delta);
        }
    }
   
};




/*  $Id: selection.js 38325 2017-04-25 22:08:59Z borodine $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 *
 * File Description: Mouse-over selection
 *
 */

Ext.define('SeqView.SelectionToolTip', {
    extend: 'NCBIGBObject.ToolTip',
    dismissDelay: 0,
    autoScroll: true,
    cls: 'SeqViewerApp',

    initComponent : function(){
        this.gview = this.selection.m_View;
        this.area = this.selection.area;
        this.pinnable = (this.area.type & SeqView.AreaFlags.NoPin) == 0;
        if (this.area.type & SeqView.AreaFlags.Track) { //track ToolTip
            this.pinnable = false;
            //this.width = 180;
            var tCfg = this.selection.m_View.m_App.m_Config.TrackConfig;
            var sig = this.area.signature;
            if (sig.indexOf(';') >= 0) sig = sig.slice(0, sig.indexOf(';'));
            for (var i = 0; i < tCfg.length; i++) {
                if (sig == tCfg[i].name) {
                    this.track = tCfg[i];
                    this.track.idx = i;
                    break;
                }
            }
        }
        this.areaShown = true;
        this.target = this.selection.element;
        this.title = this.area.title || this.selection.title;
        this.delayBeforeShow = this.showDelay;
        this.listeners = {
            beforeshow: function() {
                if (!this.delayBeforeShow) return true;
                Ext.defer(function(){
                    if (this.target) {
                        this.delayBeforeShow = 0;
                        this.show();
                    }
                }, this.showDelay, this);
                return false;
            }
        };
        this.callParent(arguments);

        if (this.track) {
            this.addTool([
                { type: 'gear',  callback: this.trackTTHandler, qtip: 'Modify Track Settings', stat: '8-1', scope: this,
                  hidden: !(this.track.check_boxes || this.track.choice_list)},
                { type: 'restore', callback: this.trackTTHandler, qtip: 'Track standalone view', stat: '8-4', scope: this },
                { type: 'help', callback: this.trackTTHandler, qtip: 'Track legend', stat: '8-2', scope: this },
                { type: 'close', hidden: this.gview.m_App.m_NoConfDlg,
                    callback: this.trackTTHandler, qtip: 'Hide (deactivate) Track', stat: '8-3', scope: this }]);
        } else {
            if (this.isPinnable()) {
                this.pinCallbackFn = function() { this.gview.pingClick('6-1'); }
                this.addTool([
                    { type: 'search', callback: this.showFeat, hidden: true, scope: this },
                    { type: 'magnify', callback: this.loadAreaRange, scope: this },
                    { type: 'up', itemId: 'up', callback: this.collapseHandler, hidden: true, scope: this },
                    { type: 'down', itemId: 'down', callback: this.collapseHandler, hidden: true, scope: this }]);
                }
        }
    },

    collapseHandler: function(o, tool, e) {
        this.collapsedTT = tool.type == 'up';
        this.setTTContent();
        this.show(this.getXY());
        //if (!this.collapsedTT) 
        this.giew.pingClick('6-4');
    },

    trackTTHandler: function(o, tool, e) {
        var app = this.selection.m_View.m_App;
        this.hide();
        this.gview.pingClick(tool.stat);
        switch (tool.type) {
            case 'close':
                this.track.shown = false;
                if (this.gview.m_SelectedSig && !this.gview.m_featMarkers)
                    SeqView.TM.Common.updateSeqViewApp(SeqView.TM.processTracksInfo(app.m_Config.TrackConfig), app);
                else 
                    this.gview.hideTrack(this.track.idx);
                break;
            case 'gear': SeqView.TM.modifyTrackDetails(this.selection.m_View, this.track); break;
            case 'help': SeqView.showHelpDlg('legends/#' + this.track.legend_text); break;
            case 'restore':
                var tracks = [this.track];
                Ext.each(app.m_Config.TrackConfig, function(trk) {
                    if (trk.shown && trk.id == "STD1") {
                        tracks.unshift(trk);
                        return false;
                    }
                });
                app.getLinkToThisPageURL('full', SeqView.TM.tracksArrayToString(tracks));
                break;
        }
    },
/*
    onShow: function() {
        this.callParent(arguments);
        if (!this.selection.area.tooltip) this.mask('Loading...');

    },
*/

    setTTContent: function() {
        var collapsedToolTips = 'NCBI/SV/ToolTips/collapse';
        if (typeof this.collapsedTT == 'undefined' && typeof localStorage !== 'undefined')
            this.collapsedTT = !(localStorage[collapsedToolTips] === undefined);
        else
            if (typeof localStorage !== undefined) {
                if (!this.collapsedTT) localStorage.removeItem(collapsedToolTips);
                else localStorage.setItem(collapsedToolTips, '');
            }

        var tt = this.area.tooltip;
        if (!tt) {
            this.mask('Loading...');
//            this.update('<div style="white-space:nowrap;">Data loading...</div>');
            this.adjustWidth = false;
            return;
        }
        if (this.isMasked() && this.getEl()) this.unmask();

        if (typeof tt == 'string') {
            this.update(tt);
            return;
        }

        if (tt.title) {this.setTitle(''); this.setTitle(tt.title)};

        var content = tt.mainInfo + (tt.downloadLinks ? tt.downloadLinks : '');

        if (tt.basicLinks) content += '<br/><b>Links & Tools</b><br/>' + tt.basicLinks;

        var tUp, tDown;
        if (tt.extraLinks) {
            tUp = this.header.down('#up');
            tDown = this.header.down('#down');
            if (tUp) tUp.hide();
            if (tDown) tDown.hide();
            if (!this.collapsedTT) {
               if (!tt.basicLinks) content += '<br/><b>Links & Tools</b>';
               content += '<br/>' + tt.extraLinks;
               tDown = false;
            } else tUp = false;
        }
        content += tt.primeButts ? tt.primeButts.html.replace(/selectionID/g, this.id) : '';

        this.update(content);
        if (tUp) tUp.show();
        if (tDown) tDown.show();
    },

    onRender: function() {
        if (!this.area.tooltip) {
            SeqView.App.simpleAjaxRequest(this.selection.ajaxCfg);
        }
        this.callParent(arguments);
    },

    onHide: function() {
        this.callParent(arguments);
        this.delayBeforeShow = this.showDelay;
        if (this.pinned) {
            delete this.selection.area.tooltip.selectionTT;

        }
        this.destroy();
    },


    afterRender : function() {
        this.callParent(arguments);
        this.setTTContent()
    },

    areaVisible: function(visible) {
        this.areaShown = visible;
    },
    showFeat: function(e,t,panel) {
        if (this.areaShown) {
            var view_start = -this.gview.m_ScrollPix;
            if (this.area.bounds.r - 10 < view_start) {
                this.gview.scrollViewTo(-(this.area.bounds.l-10), SeqView.PAN_RIGHT);
            } else  if (this.area.bounds.l+10 > this.gview.getScreenWidth() + view_start) {
                var new_pos = view_start + (this.area.bounds.r - (view_start + this.gview.getScreenWidth()));
                this.gview.scrollViewTo(-(new_pos+10), SeqView.PAN_LEFT);
            }
            this.highlightFeat();
        } else {
            this.loadAreaRange(false);
        }
        this.gview.pingClick('6-2');
    },

    highlightFeat: function() {
        var new_sel = new SeqView.SelectionHighlighter(this.gview, this.area);
        Ext.defer(new_sel.remove, 900,new_sel);
    },
    loadAreaRange: function(magnify) {
        var from = this.area.range[0];
        var len = this.area.range[1] - this.area.range[0] + 1;
        var la = Math.floor(len/10);
        from -= la;
        if (from - la < 0) from = 0;
        len += 2*la;
        if (from+len-1 > this.gview.m_App.m_SeqLen)
            len -= this.gview.m_App.m_SeqLen - (from+len-1);
        if (this.gview.m_VisFromSeq != from || this.gview.m_VisLenSeq != len) {
            this.gview.startImageLoading(from, len, {from_ui: true});
            this.delayedHighlight = true;
        } else {
            this.highlightFeat();
        }
        if (magnify) this.gview.pingClick('6-3');
    },

    highlightToolTip: function() {
        this.getEl().addCls('xsv-tooltip_highlight');
    },
    unHighlightToolTip : function(){
        this.getEl().removeCls('xsv-tooltip_highlight');
    }

});



/* Selections highlighter */

SeqView.SelectionHighlighter = function(view, area, for_edit) {
    this.m_View = view;
    this.area = area;

    this.parent_id = this.m_View.m_DivId;
    this.parent_elem = Ext.get(this.parent_id);

    var html = '<div id="selection_id_{idx}" class="sv-drag '
    html += for_edit? ' sv-highlight selection_highlight-edit' : ' selection_highlight'
    html += '"/>'
    var tpl = new Ext.Template(html);
    this.element = tpl.append(this.parent_elem, {idx: Ext.id()/*this.m_View.m_Idx*/}, true);

    var bounds = this.area.bounds;

    var off = for_edit? 6 : 3
    var left  = Math.min(view.m_Width, (view.getFlip() ? bounds.r : bounds.l) + 1);
    var the_top = bounds.t - 6;
    var the_left  = Math.max(-view.m_Width, left + this.m_View.m_ScrollPix - 4);
    var the_width = Math.min(view.m_Width * 2, Math.abs(bounds.l - bounds.r) + 6);
    var the_height = bounds.b - bounds.t + 5;

    this.element.setLeft(the_left); // the output from cgi seems a bit off. That's why this offsets. Or is it browser-specific? :)
    this.element.setTop(the_top);
    this.element.setHeight(the_height);
    this.element.setWidth(the_width);
};

SeqView.SelectionHighlighter.prototype = {
    movePix: function(delta) {
        var new_left = this.element.getLeft(true) - delta;
        this.element.setLeft(new_left);
    },
    remove: function() {
        this.element.remove();
    }
};


SeqView.createSelectionElement = function(view, area, tpl) {
    var bounds = area.bounds;
    var parent_elem = Ext.get(view.m_DivId);
    var element = tpl.append(parent_elem, {idx: Ext.id()}, true);

    var left  = Math.max(0, view.getFlip() ? bounds.r : bounds.l);
    var right  = Math.min(view.m_Width, view.getFlip() ? bounds.l : bounds.r);
    var the_width = Math.min(view.m_Width, right - left) + 2;
    var the_top = bounds.t - 1;
    var the_left  = left + view.m_ScrollPix - 1;

    if (area.type && (area.type & SeqView.AreaFlags.Track)) {
        the_left = 0;
        the_width = parent_elem.getWidth();
    }
    if (area.tname) {
        the_top = bounds.t;
        the_width = bounds.r - bounds.l;
    }
    element.setLeft(the_left); // the output from cgi seems a bit off. That's why this offsets. Or is it browser-specific? :)
    element.setTop(the_top);
    element.setHeight(bounds.b-bounds.t);
    element.setWidth(the_width);
    return element;
}
/* Selections */
SeqView.Selection = function(view, areas) {
    this.m_View = view;
    this.delayTimeout = null;
    this.saved_event = null;
    this.areas = areas;
    this.area = areas.length ? areas[areas.length - 1] : areas;

    var signatures = [];
    for (var i = 0; i != areas.length; i++) {
        signatures.push( areas[i].signature );
        if (areas[i].feats_count) {
            this.area.merged_feats = areas[i].feats_count;
        }
    }

    var tpl = (this.area.type & (SeqView.AreaFlags.NoHighlight | SeqView.AreaFlags.Track | SeqView.AreaFlags.Sequence)) !== 0 ?
            new Ext.Template('<div id="selection_id_{idx}" class="over_noselection sv-drag sv-dblclick"/>') :
            new Ext.Template('<div id="selection_id_{idx}" class="over_selection sv-drag sv-dblclick"/>');

    this.element = SeqView.createSelectionElement(view, this.area, tpl);

    // we need this complex logic for handling onClick and onDblClick events here because
    // dblclick event triggers different event sequence on different browsers:
    // on IE it is : click, dblclick
    // on others :   click, click, dblclick
    var events_handler = {
        'click':       this.clickTest,
        'contextmenu': this.onContextMenu,
        'mouseout':    this.remove,
        scope:         this
    };

//    if (Ext.isIE) {
        events_handler['dblclick'] = this.interOnDblClick;
//    }

    this.element.on(events_handler);
    var area = this.area;
    var el = this.element;
    if (!area.tooltip) {
        if (area.type & (SeqView.AreaFlags.Sequence
                       | SeqView.AreaFlags.Ruler
                       | SeqView.AreaFlags.TooltipEmbedded)) {
            area.tooltip = "Not provided";
        }
        if ((area.type & SeqView.AreaFlags.Track) !== 0 && !area.tooltip && !area.signature.match(";")){
            area.tooltip = this.m_View.m_App.getToolTipForTrack(this.m_View, area.signature);
        }
        if (area.tooltip && area.tooltip.match("Not provided")) return;
    }


    var params = null;
    if ((area.type & SeqView.AreaFlags.Track) && area.signature.match(";")){
        var gijunk = area.signature.split(";")[1];//[1] is used for getting track name
        var id = gijunk.split("|");
        params = {signatures: id[1]};
    } else {
        params = {signatures: signatures.join(',')};
    }

    var tt = area.tooltip = area.tooltip || this.m_View.m_App.findSelectionToolTip(params.signatures);
    if (!area.title) {
        if (area.range)
            area.title = (area.range[0] + 1) + ".." + (area.range[1] + 1);
        else
            area.title = area.name;
    }
    if (tt) { // use from cache
        if (tt.selectionTT && tt.selectionTT.pinned) {
           this.m_ToolTip = tt.selectionTT;
           this.m_ToolTip.selection = this;
           this.m_ToolTip.setTarget(this.element);
           this.m_ToolTip.highlightToolTip();
        }
        else this.m_ToolTip = new SeqView.SelectionToolTip({ selection: this });
        return;
    }

    params.track = this.m_View.m_App.m_Config.getObjIdTrack(area.parent_id);

    var app = this.m_View.m_App;
    if (app.m_Key) params.key = app.m_Key;
    if (app.m_DepthLimit) params.depthlimit = app.m_DepthLimit;
    if (app.m_BamPath) params.bam_path = app.m_BamPath;
    if (app.m_SRZ) params.srz = app.m_SRZ;
    if (app.m_AssmContext) params.assm_context = app.m_AssmContext;
    params.objinfo = 1;
    params.id = app.GI;
    params.color = app.m_Config.Options.curr_color;
    params.appname = this.m_View.m_App.m_AppName;
    if (area.merged_feats) {
        params.merged_feats = area.merged_feats;
    }
    this.m_ToolTip = new SeqView.SelectionToolTip({ selection: this });
    this.ajaxCfg = { url: this.m_View.m_App.m_CGIs.ObjInfo,
                    data: params, context: this,
                    success: this.checkJobStatus };
};


SeqView.Selection.prototype = {
    movePix: function(delta) {
        if (this.removed) return;
        var new_left = this.element.getLeft(true) - delta;
        this.element.setLeft(new_left);
    },

    remove: function() {
        if (this.removed) return;
        if (this.m_ToolTip) {
            if (this.m_ToolTip.pinned) this.m_ToolTip.unHighlightToolTip();
            if (this.area.tooltip) this.area.tooltip.selectionTT = this.m_ToolTip;
            this.m_ToolTip.setTarget();
        }
        this.element.remove();
        this.removed = true;
    },

    onClick: function(e) {
        if (SeqView.scrolled) {
            clearTimeout(SeqView.scrolled);
            delete SeqView.scrolled;
            return;
        }
        if (!this.m_ToolTip.pinned) this.m_ToolTip.destroy();
        this.m_View.changeSelectedSig(this.area, e.ctrlKey);
    },

    onDblClick: function(e) {
        //for ruler and track names area boundaries are ignored and zoom ins are based on coordinates
        if (this.area.type & (SeqView.AreaFlags.Track | SeqView.AreaFlags.Ruler)) return;

        var range = this.area.range;
        var new_from  = (range[0] < range[1]) ? range[0] : range[1];
        var new_to    = (range[0] < range[1]) ? range[1] : range[0];
        var new_len = new_to - new_from + 1; // before offset
        var l_offset = new_len / 100 * 6;  // calc. 5% offset
        new_from -= l_offset;
        new_to += l_offset;
        new_len = new_to - new_from + 1; // after offset
        this.m_View.startImageLoading( new_from, new_len, {from_ui: true} );
    },

    fireClick: function(e) {
        clearTimeout(this.delayTimeout);
        this.delayTimeout = null;
        if (!this.DblClickHit)
            this.onClick(this.saved_event);
        this.DblClickHit = false;
        this.saved_event = null;
    },

    clickTest: function(e) {
        if (this.delayTimeout) {
            if (e.type == 'click') {
                this.DblClickHit = false;
                clearTimeout(this.delayTimeout);
                this.delayTimeout = null;
                this.saved_event.type = 'dblclick';
                this.saved_event.browserEvent = 'Event dblclick';
                this.onDblClick(this.saved_event);
                this.saved_event = null;
            }
        } else {
            if (e.type == 'click' ){
                this.saved_event = new Ext.EventObjectImpl(e);
                var that = this;
                this.delayTimeout = setTimeout(function() {that.fireClick(e);}, 400);
            }
        }
    },

    interOnDblClick: function(e) {
        this.DblClickHit = true;
        this.onDblClick(e);
    },

    onContextMenu: function(e) {
        e.stopEvent();
        var range = this.area.range;
        if (range === undefined) return;
        var seq_pos = this.area.strand == 1 ? range[0]: range[1];
        var that = this;
        var areas;
        if (this.area.signature.indexOf('fake|') == -1) areas = this.areas;
        this.m_View.showContextMenu(e.getXY(),
            function(menu, download_menu) {
                menu.add('-', {text:'Set Sequence Origin At Feature', iconCls:'xsv-origin', scope:that,
                          handler:function() { this.m_View.m_App.showOriginDlg(seq_pos);} });
                if (areas) that.addFeatureLinks(menu, seq_pos, areas[0]);
                that.addDownloadLinks(download_menu, areas);
            }
        );
    },

    addFeatureLinks: function(menu, seq_pos, area) {
        var submenuVT = new Ext.menu.Menu();
        var menuVT = menu.add('-', {text: 'Views & Tools', iconCls: 'x-tbar-loading', disabled: true, menu: submenuVT})[1];

        var addSubmenu = function(links) {
            menuVT.setIconCls('xsv-views_tools')
            if (links.length) menuVT.enable();

            Ext.each(links, function(l) {
                submenuVT.add({text: l.label, url: l.link, handler: function(conf) { window.open(conf.url); }});
            });
        }
        if (area.links) {
            addSubmenu(area.links);
            return;
        }
        var params = this.createCGIParams(area.signature);
        params.link = 1;
        this.m_View.m_App.AjaxRequest({
            url: this.m_View.m_App.m_CGIs.Graphic, context:this, data: params,
            success: function(data, txt, res) {
                var links = SeqView.decode(data);
                for (var i = 0; i < links.length; i++) {
                    var url = links[i].link;
                    for (var j = 0; j < i; j++) {
                        if (links[j].link == url) {
                            links.splice(i--, 1);
                            break;
                        }
                    }
                }
                if (this.m_View.m_App.m_Embedded !== false)
                    for (var i = 0; i < links.length; i++) links[i].link = links[i].link.replace(/^\//, SeqView.webNCBI);
                area.links = links;
                addSubmenu(links);
            }
        }); // Ext.Ajax

    },

    addDownloadLinks: function(menu, areas) {
        if (!menu) return;
        var params = this.createCGIParams();
        var val = "";
        for (var i = 0; i < areas.length; i++) {
            if (val) val += ",";
            val += areas[i].signature;
        }
        params.ids = val;
        params.link_type = "download";
        params.link = 1;
        this.m_View.m_App.AjaxRequest({
            url: this.m_View.m_App.m_CGIs.Graphic,
            context:this, data: params,
            success: function(data, txt, res) {
                var links = SeqView.decode(data);
                var prefix = this.m_View.m_App.m_CGIs.prefix;
                for (var i = 0; i < links.length; i++) {
                    (function(url) {
                        menu.add({text:links[i].label, scope:this,
                            handler: function() {
                                var form = Ext.DomHelper.append(document.body, {
                                    tag : 'form',
                                    method : 'post',
                                    action : url
                                });
                                document.body.appendChild(form);
                                form.submit();
                                document.body.removeChild(form);
                            }
                        });
                    })(prefix + links[i].link);
                }
            }
        }); // Ext.Ajax
    },

    createCGIParams: function(signature) {
        var params = {};
        if (signature)
            params.id = signature;
        if (this.m_View.m_App.m_ViewerContext)
            params.viewer_context = this.m_View.m_App.m_ViewerContext;
        if (this.m_View.m_App.m_Key)
            params.key = this.m_View.m_App.m_Key;
        if (this.m_View.m_App.m_DepthLimit)
            params.depthlimit = this.m_View.m_App.m_DepthLimit;
        if (this.m_View.m_App.m_SRZ)
            params.srz = this.m_View.m_App.m_SRZ;
        if (this.m_View.m_App.m_BamPath)
            params.bam_path = this.m_View.m_App.m_BamPath;
        return params;
    },

    checkJobStatus: function(data, txt, req) {
        var objs = SeqView.decode(data);
        if (objs.job_status && (objs.job_status != 'failed' || objs.job_status != 'canceled')) {
            this.ajaxCfg.data = {job_key: objs.job_id};
            Ext.defer(SeqView.simpleAjaxRequest, 1000, [this.ajaxCfg]);
            return;
        }
        var obj = objs[0];
        var links = obj.links || [];
        var tip = obj.text || obj.error;
        if (objs.length > 1) { // multiple features
            tip += '<br>' + (objs[1].text || objs[1].error);
            links = links.concat(objs[1].links || []);
        }
        var app = this.m_View.m_App;
        if (app.m_preprocessorTT) {
            var prs = obj.object_id.split('-');
            app.m_preprocessorTT(links,
                {acc: app.m_Config.SeqInfo.id, title: obj.title,// GI: app.m_Config.SeqInfo.gi,
                start_pos: parseInt(prs[1], 16) + 1, end_pos: parseInt(prs[2], 16) + 1});
        }

        if (obj.exts) this.area.exts = obj.exts;
        if (obj.title) this.area.title = obj.title;

        var tt = this.area.tooltip = {mainInfo: tip, basicLinks: '', extraLinks: '', downloadLinks: '', title: this.area.title};

        var linksSet = new Set();
        for (var j = 0, jlen = links.length; j < jlen; j++) {
            var lnk = links[j];
            var path = (lnk.type == 'Download' ? SeqView.base_url : SeqView.webNCBI).slice(0, -1);
            var str = '<td align="right" valign="top" nowrap><b>';
            var link_name = (lnk.name && lnk.name.length > 0) ? lnk.name + ': ' : '';
            var link_key = window.JSON.stringify(lnk);;
            if (linksSet.has(link_key))
                continue;
            linksSet.add(link_key);
            str += link_name.replace(/ /g, '&nbsp;');
            str += '</b></td><td valign="top">';
            var anchor = lnk.html || '';
            if (!anchor) {
                for (var k = 0; k < lnk.links.length; k++) {
                    var linrec = lnk.links[k];
                    var found = false;
                    for (var l = 0; l < k; l++) {
                        if (lnk.links[l].link == linrec.link){
                            found = true;
                            break;
                        }
                    }
                    if (found) continue;

                    if (k > 0) anchor += ',&nbsp;';
                    var label = linrec.label.replace(/ /g, '&nbsp;');
                    anchor += '<a href="'
                    if (linrec.link.search(/^https?:/i) && linrec.link.search(/^ftp:/i))
                        anchor += path + (linrec.link.charAt(0) == '/' ? '' : '/') ;
                    anchor += linrec.link + '" target="_blank">' + label + '</a>';
                }
            }
            str += anchor + '</td>';
            switch (lnk.type) {
                case 'Basic': tt.basicLinks += "<tr>" + str + "</tr>"; break;
                case 'Extra': tt.extraLinks += "<tr>" + str + "</tr>"; break;
                case 'Download': tt.downloadLinks += (tt.downloadLinks ? ',&nbsp;' : '') + anchor; break;
            }
        }

        if (obj.unaligned_regions) {
            var butts = [];
            var str = '';
            var ttID = this.m_ToolTip.id;
            Ext.each(obj.unaligned_regions, function (unr, idx) {
                var id = idx + '_seq_id_' + ttID;
                butts.push({id: id, html: '<button id=\"' + id + '\" title="View unaligned fragment"'
                    + 'onclick="SeqView.clickButtonTT(\'selectionID\',' + idx + ')" '
                    + 'style=" font-size: 11px; border-radius: 4px; border-width: 1px; background-color: whitesmoke;">'
                    + 'View</button>'});

                str += '<tr><td align="right" valign="top" nowrap><b>';
                str += unescape(unr.label) + ':&nbsp;';
                str += '</b></td><td valign="top">';
                str += (unr.to - unr.from);
                str += (unr.polya ? '&nbsp;(polyA)' : '');
                str += '</td><td>&nbsp;&nbsp;&nbsp;' + butts[idx].html + '</td></tr>';
            });
            tt.primeButts = {html: '<br/><b>Unaligned Regions&nbsp;</b><table border="0" cellpadding="0" cellspacing="0"><tbody>'
                + str + '</tbody></table>', unr: obj.unaligned_regions};
        }
        if (tt.basicLinks) tt.basicLinks = '<table border="0" cellpadding="0" cellspacing="0"><tbody>' + tt.basicLinks + '</tbody></table>';
        if (tt.extraLinks) tt.extraLinks = '<table border="0" cellpadding="0" cellspacing="0"><tbody>' + tt.extraLinks + '</tbody></table>';
        if (tt.downloadLinks) tt.downloadLinks = '<br/><b>Download:&nbsp;</b>' + tt.downloadLinks + '<br/>';

        this.m_View.m_App.saveSelectionToolTip(this.ajaxCfg.data.signatures, tt);

        if (this.m_ToolTip.isVisible()) this.m_ToolTip.setTTContent();

        this.m_View.m_App.ping({"jsevent": "mouseover","sv-event":"Object_ToolTip"});
    },

    showPrime: function(unr, params) {
        var urHead = (unr.label[0] == '5');
        var reverse = unr.orientation == 'reverse';
        var length = Math.abs(unr.to - unr.from);
        var tail = Math.min(unr.length - (unr.to - unr.from), 500);
        if (length > 2000) {
            if (urHead + reverse == 1) unr.from = unr.to - 2000;
            else unr.to = unr.from + 2000;
        }
        params = params || {url: this.m_View.m_App.m_CGIs.Graphic, context: this,
            data: {objinfo: 1, id: unr.seq_id, unalignedregion: true, from: unr.from, to: unr.to,
                alignedtail: tail, reverse: reverse}};
        if (this.m_View.m_Flip != params.data.reverse) urHead = ! urHead;
        params.success = function(data){
            switch (data.job_status) {
                case 'submitted': case 'running': case 'pending':
                    params.data.job_id = data.job_id;
                    Ext.defer(this.showPrime, 1000, this, [unr, params]);
                    break;
                default:
                    var up = data.unaligned_part,
                        ap = data.aligned_part;
                    var rgText = '<p style="word-wrap:break-word; margin: 3px;">' + (urHead ? '' : ap.sequence) + '<mark style="background-color:red">' + up.sequence + '</mark>' + (urHead ? ap.sequence : '') + '</p>';
                    var mbw = new Ext.Window({
                        title: 'Unaligned Region (' + unr.orientation + ' strand)',
                        app: this.m_View.m_App,
                        modal: true, layout:'fit',
                        overflowY: 'auto',
                        width: 400, height: 300,
                        items:[{xtype: 'displayfield', value: rgText, textalign: 'left'}],
                        buttons:[{text: 'Close', handler: function() {mbw.close();}}]});
                    mbw.show();
            }
        }
        params.error = function(d, t){
            Ext.MessageBox.alert('Error', d);
        }
        this.m_View.m_App.AjaxRequest(params);
    }
};

SeqView.clickButtonTT = function(id, idx) {
    var selTT = Ext.getCmp(id);
    if (selTT) selTT.selection.showPrime(selTT.area.tooltip.primeButts.unr[idx]);
};



Ext.define('SeqView.SelectedRangeToolTip', {
    extend: 'NCBIGBObject.ToolTip',
    anchor: 'top',
    anchorToTarget: true,
    initComponent: function() {
        this.gview = this.selection.m_View;
        this.target = this.selection.m_tt_TargetDivId;
        this.id = this.selection.m_tt_DivId;
        this.title = 'Range: ' + this.selection.getRangeStr();

        this.callParent(arguments);
        var mHeight = 0;
        var nextStep = function () {
            return mHeight += 22;
        };
        var toolbar_items = [{
            text: 'Zoom On Range', iconCls: 'xsv-zoom_range',
            style: {margin: '0px',padding:'0px'}, x: 0, y: mHeight, scope: this,
            handler: function() {
                this.hide();
                var range = this.selection.range;
                this.gview.removeRangeSelection(range);
                this.gview.startImageLoading(range[0], range[1] - range[0] + 1, {from_ui: true});
                this.gview.pingClick('5-0');
            }},{
            text:'Zoom To Sequence', iconCls: 'xsv-zoom_seq',
            style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
            handler:function() {
                this.hide();
                var range = this.selection.range;
                this.gview.removeRangeSelection(range);
                this.gview.zoomSeq(Math.floor((range[1] + range[0])/2));
                this.gview.pingClick('5-8');
            }},{
            text: 'Modify Range', iconCls: 'xsv-zoom_range',
            style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
            handler: function() {
                this.hide();
                this.modifyRange();
                this.gview.pingClick('5-1');
            }}];

        if (this.gview.m_App.mf_MultiPanel) {
            toolbar_items.push({
                text: 'Add New Panel On Range', iconCls: 'xsv-new_view',
                style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
                handler: function() {
                    var view = new SeqView.Graphic(this.gview.m_App);
                    this.gview.m_App.registerView(view);

                    var from = -1, to = -1;
                    var range = this.gview.getTotalSelectedRange();
                    if (range[0] !== -1 && range[1] !== -1) {
                        from = this.selection.range[0];
                        to   = this.selection.range[1];
                    } else if (this.gview.m_UrlFrom) {
                        from = this.gview.m_UrlFrom;
                        to   = this.gview.m_UrlTo;
                    } else {
                        from = this.gview.m_VisFromSeq + 1;
                        to   = this.gview.m_VisFromSeq + this.gview.m_VisLenSeq;
                    }

                    if (from !== -1 && to !== -1) view.startImageLoading(from, to - from + 1);
                    this.hide();
                    this.highlightTargetRectangle();
                    this.gview.pingClick('5-2');
                }}
            );
        }
        var noPBlast = (this.gview.m_App.m_ViewParams['acc_type'] !== 'DNA');
        toolbar_items.push({
            text: 'Set New Marker For Selection', iconCls: 'xsv-markers',
            style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
            handler: function() {
                this.hide();
                this.gview.m_App.addMarker(this.selection.range);
                this.gview.removeRangeSelection(this.selection.range);
                this.gview.pingClick('5-3');
            }},{
            text: 'BLAST Search (Selection)', iconCls: 'xsv-blast',
            style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
            handler: function() {
                this.gview.blastSelection();
                this.hide();
                this.gview.pingClick('5-4');
            }}, {
            text: 'Primer BLAST (Selection)', iconCls: 'xsv-primer', disabled: noPBlast,
            style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
            handler: function() {
                this.gview.primerBlast();
                this.hide();
                this.gview.pingClick('5-5');
            }},{
            text: 'Download FASTA (Selection)', iconCls: 'xsv-download-static',
            style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
            handler: function() {
                this.gview.downloadData(true, "fasta", this.selection.range);
                this.gview.pingClick('5-6');
            }},{
            text: 'Download GenBank Flat File (Selection)', iconCls: 'xsv-download-static',
            style: {margin: '0px',padding:'0px'}, x: 0, y: nextStep(), scope: this,
            handler: function() {
                this.gview.downloadData(true, "flat", this.selection.range);
                this.gview.pingClick('5-7');
            }}
        );


        var toolbar = new Ext.Toolbar({
            layout: 'absolute',
            border: false,
            items: toolbar_items,
            height: mHeight + 16,
            width: 218
        });

        this.add(toolbar);
    },

    modifyRange: function() {
        Ext.MessageBox.prompt('Modify Range', 'Enter New Range:', function(btn, text)
        {
            if (btn!='ok' || text.length==0) return;
            this.gview.m_App.handlePos(text, {
                allow_equal: true,
                ask_user: true,
                success: function(pos_range, options) {
                    if (pos_range.length !== 2) {
                        options.failure.call(options.scope, "Range required", options);
                        return;
                    }
                    this.selection.range[0] = pos_range[0];
                    this.selection.range[1] = pos_range[1];
                    this.selection.updateCoords();
                    this.setTitle('Range: ' + this.selection.getRangeStr());
                },
                failure: function(err_msg, options) {
                    Ext.MessageBox.alert('Error', err_msg);
                },
                scope: this
            });
        }, this, false, this.selection.getRangeStr());
    },

    onTargetOut: function(e){
        if (this.disabled || e.within(this.target.dom, true)){
            return;
        }
        this.highlightTargetRectangle();
        this.el.applyStyles('outline-style: none');
    },

    onTargetOver: function(e){
        this.highlightTargetRectangleMouseOver();
        if (this.pinned) {
            this.el.applyStyles('outline-style: solid;outline-color: red;outline-width: thin;');
        }
//        else this.show();
        this.show();
//        this.gview.ping({"jsevent":"mouseover","sv-event":"Range_ToolTip"});
    },

    onShow: function() {
        if (this.gview.m_ContextMenu && this.gview.m_ContextMenu.isVisible()){
                this.gview.m_ContextMenu.hide();
         }
         this.callParent(arguments);
    },

    highlightTargetRectangleMouseOver: function() { this.target.dom.style.backgroundColor = '#254117'; },
    highlightTargetRectangle: function() { this.target.dom.style.backgroundColor = '#95B9C7'; }
});

//////////////////////////////////////////////////////////////////////////
SeqView.RangeSelection = function(view, range, event) {
    this.m_View = view;
    this.m_Resizing = false;
    if (!event) this.range = range;

    var parent_elem = Ext.fly(view.m_DivId);
    this.m_tt_TargetDivId = 'svtt-ruler-' + Ext.id();
    this.m_tt_DivId = 'svtt-' + Ext.id();
    var TargetDiv = '<div id=' + this.m_tt_TargetDivId + ' style="height:' + view.m_topRulerHeight + 'px; border: none; background-color: #95B9C7;"></div>';
    if (view.m_bottomRulerHeight)
        TargetDiv +='<div id=' + ('svtt-ruler-' + Ext.id()) + ' style="height:' + view.m_bottomRulerHeight + 'px; top:'
            + (view.m_Height - view.m_bottomRulerHeight - view.m_topRulerHeight)
            + 'px; position:relative; border:none; background-color: #95B9C7;"></div>';
    var RangeDiv = '<div class="range_selection sv-drag sv-highlight" style="z-index:1;">'+ TargetDiv +'</div>';

    var tpl = new Ext.Template(RangeDiv);
    this.element = tpl.append(parent_elem, {}, true);
    this.m_ToolTip = null;

    this.element.on({
        'mousedown': this.onMouseDown,
//        'mousemove': this.onMouseMove,
        'touchdown': this.onMouseDown,
        'touchmove': this.onMouseMove,
        'click':     this.onClick,
        scope: this
    });
    if (view.m_bottomRulerHeight) this.element.on({'mouseover': this.onMouseOver, scope: this});

     if (event /*means resizing*/) {
        var x = this.pageToViewX(Math.round(range[0]));
        this.range = [x, x];
        this.mouse_offset = 0;
        this.setElementPixRange(this.range);
        this.startResizing(event);
    } else {
        this.updateCoords();
    }
    if (!event && !this.m_ToolTip) {
        this.m_ToolTip = new SeqView.SelectedRangeToolTip({selection: this});
    }
};

SeqView.RangeSelection.prototype = {
    pageToViewX: function(x) {
        var div_xy = Ext.fly(this.m_View.m_DivId).getXY();
        return x - div_xy[0];
    },

    getRangeStr: function() {
        var app = this.m_View.m_App;
        var start = app.posToLocalDisplay(this.range[0]);
        var end = app.posToLocalDisplay(this.range[1]);
        return start + '..' + end;
    },

    movePix: function(delta) {
        var new_left = this.element.getLeft(true) - delta;
        this.element.setLeft(new_left);
    },

    remove: function() {
        this.element.remove();
        if (this.m_ToolTip) this.m_ToolTip.hide();
    },

    rangePixToSeq: function(range) {
        var view = this.m_View;
        var scp = view.m_ScrollPix;
        var over = view.getFlip() ? -0.5 : 0.5;
        var left = Math.round(view.pix2Seq(range[0]-scp, over));
        var right = Math.round(view.pix2Seq(range[1]-scp, -over));
        if (left <= right) 
            return [left, right];
        else
            return [right, left];
    },

    rangeSeqToPix: function(range) {
        var view = this.m_View;
        var left, right, over;
        if (!view.getFlip()) {
            over = 0.5;
            left = range[0];
            right = range[1];
        } else {
            over = - 0.5;
            left = range[1];
            right = range[0];
        }
        left = view.seq2PixScrolled(left - over);
        right = view.seq2PixScrolled(right + over);
        return [left, right];
    },

    updateCoords: function() {
        var pix_range = this.rangeSeqToPix(this.range);
        this.setElementPixRange(pix_range);
    },
    setElementPixRange: function(range) {
        var view = this.m_View;
        // For logic underneath see markers.js:updateMarkerPos
        // TODO: move this code to view
        var left = Math.max(-10000, Math.min(10000, range[0]));
        var right = Math.max(-10000, Math.min(10000, range[1]));
        var width = right - left;
        if (width > 4000) {
            var mid_view = Math.round(view.getScreenWidth() / 2);
            if (Math.abs(left - mid_view) < Math.abs(right - mid_view)) {
                right = Math.max(mid_view + 2000, Math.min(right, left + 4000));
                left = Math.max(left, right - 4000)
            } else {
                left = Math.min(mid_view - 2000, Math.max(left, right - 4000));
                right = Math.min(right, left + 4000);
            }
            width = right - left;
        }
        var el = this.element;
        el.setLeft(left);
        el.setWidth(width);
        el.setTop(view.getSelectionTop());
        el.setHeight(view.getSelectionHeight());
    },

    startResizing: function(e) {
        this.m_Resizing = true;
        var o = this;
        var onMove = function(e) {o.onMouseMove(new Ext.event.Event(e ? e : window.event));};
        var onEnd = function(e) {o.onMouseUp(new Ext.event.Event(e ? e : window.event));};
        if ('button' in e && e.button == 0) {
            this.m_DocMouseMove = document.onmousemove;
            this.m_DocMouseUp = document.onmouseup;
            document.onmousemove = onMove;
            document.onmouseup = onEnd;
        } else {
            this.m_DocTouchMove = document.ontouchmove;
            this.m_DocTouchUp = document.ontouchend;
            document.ontouchmove = onMove;
            document.ontouchend = onEnd;
        }
    },

    onClick: function(e) {
        if (e.ctrlKey) {
            e.stopEvent();
            this.m_View.removeRangeSelection(true, this);
        }
    },

    onMouseOver: function(e) {
        var xy = e.getXY();
        if (this.m_ToolTip && e.target.id.indexOf('svtt-ruler') == 0) {
            if (this.m_ToolTip.target.id != e.target.id) {
                this.m_ToolTip.setTarget(e.target.id);
                this.m_ToolTip.show();
            }
        }
    },

    onMouseDown: function(e) {
        // if hit is close to the border, turn on resize
        var xy = e.getXY();
        var x = this.pageToViewX(xy[0]);
        var left = this.element.getLeft(true);
        var right = left + this.element.getWidth(true);
        var ld = Math.abs(x - left);
        var rd = Math.abs(x - right);
        if (ld < 5 || rd < 5) {
            e.stopPropagation();
            var pix_range = this.rangeSeqToPix(this.range);
            if (ld < rd) this.range = [pix_range[1], pix_range[0]];
            else this.range = pix_range;
            this.mouse_offset = x - this.range[1];
            this.m_View.popRangeElement(this);
            this.startResizing(e);
        }
    },

    onMouseUp: function(e) {
        var range;
        if (this.range[0] > this.range[1])
            range = [this.range[1], this.range[0]];
        else
            range = this.range;
        
        this.range = this.rangePixToSeq([Math.max(0, range[0]), Math.min(range[1], this.m_View.m_Width)]);
        this.updateCoords();
        if (this.m_Resizing) {
            Ext.fly(this.m_View.m_DivId).setStyle('cursor', 'default');
            Ext.fly(e.getTarget()).setStyle('cursor', 'default');
        }

        this.m_Resizing = false;

        if ('button' in e && e.button == 0) {
            document.onmousemove = this.m_DocMouseMove;
            document.onmouseup = this.m_DocMouseUp;
        } else {
            document.ontouchmove = this.m_DocTouchMove;
            document.ontouchend = this.m_DocTouchUp;
        }
        this.m_DocMouseMove = this.m_DocMouseUp = this.m_DocTouchMove = this.m_DocTouchUp = null;

        // Let view decide whether it needs this range selection
        if (!this.m_View.addRangeSelection(this, this.range, range[1]-range[0])) return;

        if (!this.m_ToolTip)
            this.m_ToolTip = new SeqView.SelectedRangeToolTip({ selection: this, autoHide: true });
        if (this.m_ToolTip) {
            var el = Ext.get(this.m_tt_TargetDivId);
            if (el) {
                this.m_ToolTip.show();
                this.m_ToolTip.setTitle('Range: ' + this.getRangeStr());
                this.m_ToolTip.highlightTargetRectangleMouseOver();
            }
            this.m_View.m_App.ping({"jsevent":"click","sv-event":"Graph_Range_Selection"});
        }
    },

    onMouseMove: function(e) {
        NCBIGBUtils.ClearBrowserSelection();
        var xy = e.getXY();
        var x = this.pageToViewX(xy[0]);
        if (!this.m_Resizing) {
            var left = this.element.getLeft(true);
            var right = left + this.element.getWidth(true);
            var ld = Math.abs(x - left);
            var rd = Math.abs(x - right);
            if (ld < 5 || rd < 5) {
                Ext.fly(e.getTarget()).setStyle('cursor', 'e-resize');
            } else {
                Ext.fly(e.getTarget()).setStyle('cursor', 'default');
            }
            if (this.m_ToolTip && !this.m_ToolTip.isVisible()) {
                this.m_ToolTip.highlightTargetRectangle();
            }
            return;
        }
        if (e.type == 'mousemove' && e.button < 0) { this.onMouseUp(e);
        } else {
            x -= this.mouse_offset;
            e.stopEvent();
            if (x == this.range[1]) return;
            this.range[1] = x;
            var left, width;
            if (this.range[0] < this.range[1]) {
                this.setElementPixRange(this.range);
            } else {
                this.setElementPixRange([this.range[1], this.range[0]]);
            }
        }
    }
};/*  $Id: markers.js 38325 2017-04-25 22:08:59Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 * File Description: Visual markers
 *
 */

SeqView.MarkerFlags = {
    UserLock:      (1 << 0), ///< Locked by user
    SystemLock:    (1 << 1), ///< Locked by Sequence Viewer, unlockable by user
    Hollow:        (1 << 2), ///< Marker body is transparent
    PositionTitle: (1 << 3)  ///< Title is position-based and cleared if pos change
};
SeqView.MarkerNav = String.fromCharCode(171, 187);

//////////////////////////////////////////////////////////////////////
// SeqView.Marker

SeqView.Marker = (function () {
    var sm_Width = 8;

    // Coordinates are zero-based sequence position
    function constructor(minfo, range, name, flags, color) {
        this.m_MInfo = minfo;
        this.obj_coords = null;
        this.deleted = false;
        this.flags = flags;

        this.color = color || minfo.getNextColor();
        this.assignName(name ? name : minfo.suggestNextName());
        this.marker_num = minfo.m_NextMarkerNum++;

        this.span = 0;
        this.seq_pos = range[0];
        if (range[1])
            this.span = range[1] - range[0] + 1;

        this.marker_in_views = [];
        this.lines_in_views = [];
        this.addToAllViews();


    };

    constructor.getWidth = function() { return sm_Width; };

    var sm_PanoramaTmpl = new Ext.Template(
      '<div id="{id}" data-qtip="{qtip}" class="sv-marker sv-marker_base sv-opaque67">',
      '<div style="position:absolute;top:7px;left:0px;width:16px;height:2px;clip:rect(0 16px 2px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:4px;left:1px;width:14px;height:8px;clip:rect(0 14px 8px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:3px;left:2px;width:12px;height:10px;clip:rect(0 12px 10px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:2px;left:3px;width:10px;height:12px;clip:rect(0 10px 12px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:1px;left:4px;width:8px;height:14px;clip:rect(0 8px 14px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:0px;left:7px;width:2px;height:16px;clip:rect(0 2px 16px 0);background-color:{color};"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:1px;border-left:1px solid {color};">',
          '<div id="{id}_body" class="sv-opaque3 sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div data-qtip="Locked" id="the_{id}_lock" style="display:{lock};" class="marker_locked_ov"></div></div>'
    );
    
    var sm_PanoramaRectTmpl = new Ext.Template(
      '<div id="{id}" data-qtip="{qtip}" class="sv-marker sv-marker_base sv-opaque67">',
      '<div style="position:absolute;top:-2px;left:2px;width:13px;height:14px;background-color:{color};opacity:0.4;"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="top:14px;width:1px;border-left:1px solid {color};">',
      '</div></div>'
    );
    

    var sm_GraphicTmpl = new Ext.Template(
      '<div id="{id}" class="sv-marker sv-marker_half_base" style="top:0px;">',
      '<div style="position:absolute;top:7px;left:0px;width:8px;height:2px;clip:rect(0 8px 2px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:4px;left:1px;width:7px;height:8px;clip:rect(0 7px 8px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:3px;left:2px;width:6px;height:10px;clip:rect(0 6px 10px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:2px;left:3px;width:5px;height:12px;clip:rect(0 5px 12px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:1px;left:4px;width:4px;height:14px;clip:rect(0 4px 14px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:0px;left:7px;width:1px;height:16px;clip:rect(0 1px 16px 0);background-color:{color};"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:0px;border-left:1px solid {color};border-right:1px solid {color};">',
          '<div id="{id}_body" class="sv-opaque3 sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div id="the_{id}_back" class="marker_label_back sv-opaque67" style="background-color:{color};width:{label_length}px;"></div>',
       '<div id="{id}_label" class="marker_label">{trimmed_label}',
      ' <img id="the_{id}_lock" style="margin-top:1px;display:{lock};" src="' + SeqView.base_url + 'images/lock.png"></div></div>'
    );
    
    var sm_GraphicRectTmpl = new Ext.Template(
      '<div id="{id}" class="sv-marker sv-marker_half_base" style="top:0px;">',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:0px;border-left:1px solid {color};border-right:1px solid {color};">',
          '<div id="{id}_body" class="sv-opaque3 sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div id="the_{id}_back" class="marker_label_back sv-opaque67" style="background-color:{color};opacity:0.4;top:1px;left:1px;width:18px;"></div>',
      '<div id="{id}_label" class="marker_label" style="font-size:13px;left:1px;top:0px;">' + SeqView.MarkerNav + '</div>',
      '</div>'
    );


    var sm_GraphicHollowTmpl = new Ext.Template(
      '<div id="{id}" class="sv-marker sv-marker_half_base" style="top:0px;">',
      '<div style="position:absolute;top:7px;left:0px;width:8px;height:2px;clip:rect(0 8px 2px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:4px;left:1px;width:7px;height:8px;clip:rect(0 7px 8px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:3px;left:2px;width:6px;height:10px;clip:rect(0 6px 10px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:2px;left:3px;width:5px;height:12px;clip:rect(0 5px 12px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:1px;left:4px;width:4px;height:14px;clip:rect(0 4px 14px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:0px;left:7px;width:1px;height:16px;clip:rect(0 1px 16px 0);background-color:{color};"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:0px;border-top:5px solid {color};border-left:1px solid {color};border-right:1px solid {color};">',
          '<div id="{id}_body" class="sv-transparent sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div id="the_{id}_back" class="marker_label_back sv-opaque67" style="background-color:{color};width:{label_length}px;"></div>',
      '<div id="{id}_label" class="marker_label">{trimmed_label}',
      ' <img id="the_{id}_lock" style="margin-top:1px;display:{lock};" src="' + SeqView.base_url + 'images/lock.png"></div></div>'
    );

    constructor.getTemplate = function(is_panorama, fill_bg, rect) {
        if (is_panorama) return rect ? sm_PanoramaRectTmpl : sm_PanoramaTmpl;
        if (rect) return sm_GraphicRectTmpl;
        return fill_bg ? sm_GraphicTmpl : sm_GraphicHollowTmpl;
    };

    return constructor;
}) ();


SeqView.Marker.prototype = {

    displayPos: function(forPDF) {
        var pos = forPDF ? this.seq_pos : this.seq_pos + 1;
        if (this.span) pos += ':' + (this.seq_pos + this.span);
        return pos;
    },

    displayLocalPos: function() {
        var app = this.m_MInfo.m_App;
        var pos = app.posToLocal(this.seq_pos);
        if (this.span) {
            if (app.getFlip()) {
                pos = app.posToLocalDisplay(this.seq_pos + this.span - 1) + ':' + pos;
            } else {
                pos += ':' + app.posToLocalDisplay(this.seq_pos + this.span - 1);
            }
        }
        return pos;
    },

    addToAllViews: function() {
        var app = this.m_MInfo.m_App;
        app.forEachView(function(view) { this.addToView(view); }, this);
        if (app.m_TextView)
            app.m_TextView.AddOrUpdateMarker(this); // text view
        app.fireEvent('marker_created', this);
    },

    getCreateParams: function(is_panorama, idx) {
        var elem_id = 'marker_' + idx + '_' + this.marker_num;
        var user_lock = this.flags & SeqView.MarkerFlags.UserLock;
        this.nvMarker = this.marker_name == SeqView.MarkerNav;
        return {
            template: SeqView.Marker.getTemplate(is_panorama,
                !(this.flags & SeqView.MarkerFlags.Hollow), this.nvMarker),
            options: {
                id: elem_id, color: this.color,
                num: this.marker_num, trimmed_label: this.marker_name_trimmed,
                label_length: user_lock ? this.marker_name_back_lock : this.marker_name_back,
//                font_size: this.m_MInfo.m_App.m_Portal ? '0.85em' : '11px',
                lock: user_lock ? 'inline' : 'none', prefix: this.m_MInfo.m_App.m_CGIs.html_prefix
            }
        };
    },

    addToView: function(view) {
        if (view.isAlignment())
            return;

        var elem_id = 'marker_' + view.m_Idx + '_' + this.marker_num;

        var marker = view.createMarkerElem(this);
        this.m_MInfo.m_AllMarkers[elem_id] = this;

        var marker_line = Ext.get(elem_id + '_line');
        var marker_height = view.m_FromCgi.img_height || view.m_Height;
        marker_line.setHeight(marker_height - 18); // ball size
        // To satisfy IE pre-8
        marker.setHeight(marker_height);

        this.lines_in_views[view.m_Idx] = marker_line;
        this.marker_in_views[view.m_Idx] = marker;

        this.updateMarkerPos(view);

        var marker_label = Ext.get(elem_id + '_label');
        if (marker_label) {
            marker_label.on({
                'mousedown': this.onMouseDown,
                'mousemove': this.onMouseMove,
                'contextmenu': this.onContextMenu,
                scope: this
            });
        }
        //creating tooltip object that targets this marker
        if (view.m_Type == "graphical") {
            this.m_ToolTip = new SeqView.MarkerToolTip({
                selection: this,
                target: marker_label,
                gview: view
            });
        }

    },
    removeFromView: function(view) {
        //removing marker tooltips
        if (this.m_ToolTip) {
//            this.m_ToolTip.remove();
            this.m_ToolTip.destroy();
            delete this.m_ToolTip;
        }
        var elem_id = 'marker_' + view.m_Idx + '_' + this.marker_num;
        if (this.marker_in_views[view.m_Idx]) {
            this.marker_in_views[view.m_Idx].remove();
            delete this.m_MInfo.m_AllMarkers[elem_id];
        }
    },

    lockMarker: function() {
        if (this.flags & SeqView.MarkerFlags.SystemLock) return;
        this.flags = this.flags ^ SeqView.MarkerFlags.UserLock;
        if (!(this.flags & SeqView.MarkerFlags.UserLock) &&
             (this.flags & SeqView.MarkerFlags.PositionTitle)) {
            this.flags = this.flags ^ SeqView.MarkerFlags.PositionTitle;
            this.rename(this.m_MInfo.suggestNextName());
        }

        // update views
        this.m_MInfo.m_App.forEachView(function(view) {
            var elem_id = 'the_marker_' + view.m_Idx + '_' + this.marker_num + '_lock';
            var back_id = 'the_marker_' + view.m_Idx + '_' + this.marker_num + '_back';
            var the_marker = Ext.get(elem_id);
            var the_back = Ext.get(back_id);

            // it's a bit tricky. Marker backgroung is transparent and the label is opaque
            // this means we have to use absolute positioning and resize the background manually to put/remove the lock icon
            if (the_marker) {
                var cur_display = the_marker.getStyle('display');
                the_marker.setStyle('display', cur_display == 'none' ? 'inline' : 'none');
                if (!view.isPanorama()) {
                    the_back.setStyle('width', '' + ((this.flags & SeqView.MarkerFlags.UserLock) ? this.marker_name_back_lock : this.marker_name_back) + 'px');
                }
            }
        }, this);

        // Text View
        if (this.m_MInfo.m_App.m_TextView)
            this.m_MInfo.m_App.m_TextView.UpdateLockMarker(this);
    },

    deleteMarker: function() {
        this.m_MInfo.m_App.forEachView(function(view) {
            this.removeFromView(view);
        }, this);

        // remove from text view
        if (this.m_MInfo.m_App.m_TextView)
            this.m_MInfo.m_App.m_TextView.RemoveMarker(this);

        this.deleted = true;
        this.m_MInfo.m_App.fireEvent('marker_deleted', this);

    },


    getMarker: function(elem_id) {
        var s = elem_id.split('_');
        view = this.m_MInfo.m_App.findView(s[1]);
        marker = this.marker_in_views[s[1]];
        return { view: view, marker: marker };
    },

    setCursor: function(elem_id, cursor) {
        this.getMarker(elem_id).marker.setStyle('cursor', cursor);
    },

    getSeqPos: function(elem_id, pix_pos) {
        var s = elem_id.split('_');
        var view = this.m_MInfo.m_App.findView(s[1]);
        if (view) {
            var pix = pix_pos - view.m_ScrollPix;
            return view.pix2Seq(pix);
        }
        return 0;
    },

    // set marker name and calculate derived names and parameters
    assignName: function(name) {
        this.marker_name_trimmed = name.trimToPix(100);
        this.marker_name = name;
        this.visible_name = NCBIGBUtils.sanitize(name);
        this.marker_name_back = this.marker_name_trimmed.visualLength() + 2;
        this.marker_name_back_lock = this.marker_name_back + 16;
    },

    setName: function(name) {
        this.assignName(name);

        this.m_MInfo.m_App.forEachView(function(view) {
            if (view.isAlignment())
                return;

            var elem_id = 'marker_' + view.m_Idx + '_' + this.marker_num;
            var back_id = 'the_marker_' + view.m_Idx + '_' + this.marker_num + '_back';
            var elem_label_id = 'marker_' + view.m_Idx + '_' + this.marker_num + '_label';

            var marker_div = Ext.get(elem_id);
            var the_back = Ext.get(back_id);
            var marker_label_div = Ext.get(elem_label_id);
            if (marker_label_div) {// no labels in overview
                var lock_html = marker_label_div.dom.innerHTML;
                lock_html = lock_html.slice(lock_html.indexOf('<'), lock_html.length);

                marker_label_div.update(this.marker_name_trimmed + ' ' + lock_html);
                the_back.setStyle('width', '' + ((this.flags & SeqView.MarkerFlags.UserLock) ? this.marker_name_back_lock : this.marker_name_back) + 'px');
            }
        }, this);
        var marker_seq = Ext.get('seqtext-marker_' + this.m_MInfo.m_App.m_Idx + '_' + this.marker_num);
    },

    rename: function(name) {
        this.setName(name);
        //changing tooltip title too
        var title = '&nbsp;&nbsp;&nbsp;' + name;
        this.m_ToolTip.setTitle(title);
        //
        this.m_MInfo.updateInfo();
    },

    setSeqPos: function(pos, reload) {
        this.seq_pos = pos;

        this.m_MInfo.m_App.forEachView(function(view) {
            if (view.isAlignment())
                return;
            this.updateMarkerPos(view);
        }, this);

        if (reload) {
            this.m_MInfo.updateInfo();
        }

        if (reload && this.m_MInfo.m_App.m_TextView) { // text view
            this.m_MInfo.m_App.m_TextView.SetSeqPos(this);
        }
    },

    setSeqRange: function(range, reload) {
        var pos = range[0];
        if (range[1] && range[1] > pos)
            this.span = range[1] - pos + 1; // inclusive range
        else
            this.span = 0;
        this.setSeqPos(pos, reload);
    },

    setTop: function(h) { this.marker.setTop(h); },

    initMarkerMove: function(marker_id, xy) {
        this.m_CurMarker = marker_id;
        this.m_PrevXY = xy;
        var vm = this.getMarker(marker_id);
        this.m_mouseOffset = vm.marker.getLeft(true) - xy[0];
        this.setCursor(this.m_CurMarker, 'move');
    },

    moveMarker: function(xy) {
        var vm = this.getMarker(this.m_CurMarker);
        var new_left = xy[0] + this.m_mouseOffset;
        var marker_width = SeqView.Marker.getWidth();
        var flip = vm.view.getFlip && vm.view.getFlip(); // ugly, but overview does not have a getFlip
        // validate coordinate against visible window
        new_left = Math.max(-marker_width, new_left);
        var rightmost_pix = vm.view.getMostRightPix();
        var dWidth = vm.marker.getWidth() - marker_width;
        new_left = Math.min(rightmost_pix - dWidth, new_left);
        var pix = new_left + marker_width - vm.view.m_ScrollPix + ((dWidth && !this.span) ? vm.view.m_bpPix/2 : 0);
        var new_seq_pos = vm.view.pix2Seq(pix) - (flip ? this.span : 0);
        // Check new position in sequence space and fix it if possible:
        // taking into account vm.marker.getWidth() (above) takes care of this
        // approximately. Should work worse for longer sequences where
        // ration between sequence space and pixel space is higher.
        // So we adjust new_seq_pos so that both ends of marker are in
        // valid range.
        var fixed_seq_pos = Math.max(0, new_seq_pos);
        fixed_seq_pos = Math.min(vm.view.m_App.m_SeqLength - this.span, fixed_seq_pos);
        this.seq_pos = fixed_seq_pos;
        vm.marker.setLeft(new_left);
        // syncronize other markers
        this.m_MInfo.m_App.forEachView(function(view) {
            if (view.isAlignment() || view.m_Idx == vm.view.m_Idx)
                return;
            this.updateMarkerPos(view);
        }, this);

        if (this.m_MInfo.m_App.m_TextView) { // text view
            this.m_MInfo.m_App.m_TextView.AddOrUpdateMarker(this);
        }
    },

    movePix: function(elem_id, delta, validate) {
        var vm = this.getMarker(elem_id);
        var new_left = vm.marker.getLeft(true) - delta;
        var marker_width = SeqView.Marker.getWidth();
        if (validate) { // validate only when dragging with mouse
            new_left = Math.max(-marker_width, new_left);
            var rightmost_pix = vm.view.getMostRightPix();
            new_left = Math.min(rightmost_pix /* - vm.marker.getWidth() */ + marker_width, new_left);
        }
        var pix = new_left + marker_width - vm.view.m_ScrollPix;
        var flip = vm.view.getFlip && vm.view.getFlip(); // ugly, but overview does not have a getFlip
        var new_seq_pos = vm.view.pix2Seq(pix) - (flip ? this.span : 0);
        // Check new position in sequence space and fix it if possible
        var fixed_new_seq_pos = Math.max(0, new_seq_pos);
        fixed_new_seq_pos = Math.min(vm.view.m_App.m_SeqLength - this.span - 1, fixed_new_seq_pos);
        this.seq_pos = fixed_new_seq_pos;
        if (fixed_new_seq_pos !== new_seq_pos) {
            new_left = vm.view.seq2Pix(fixed_new_seq_pos);
        }
        vm.marker.setLeft(new_left);

        if (validate) { // syncronize other markers
            this.m_MInfo.m_App.forEachView(function(view) {
                if (view.isAlignment() || view.m_Idx == vm.view.m_Idx)
                    return;
                this.updateMarkerPos(view);
            }, this);

            if (this.m_MInfo.m_App.m_TextView) { // text view
                this.m_MInfo.m_App.m_TextView.AddOrUpdateMarker(this);
            }
        }
    },

    scrollPix: function(view, delta) {
        if (view.isAlignment())
            return;

        var marker = this.marker_in_views[view.m_Idx];
        var new_left = marker.getLeft(true) - delta;
        marker.setLeft(new_left);
    },

    updateMarkerSize: function(view) {
        if (view.isAlignment())
            return;
        var idx = view.m_Idx;
        var marker_line;
        if (this.lines_in_views[idx] != null) {
            marker_line = this.lines_in_views[idx];
        } else {
            this.addToView(view);
            marker_line = this.lines_in_views[idx];
        }
        marker_line.setHeight((view.m_FromCgi.img_height || view.m_Height) - 18);
    },

    updateMarkerPos: function(view) {
        if (view.isAlignment())
            return;
        var view_idx = view.m_Idx;
        var marker = this.marker_in_views[view_idx];
        var effective_pos = this.seq_pos - 0.5;
        var pix_pos = view.seq2PixScrolled(effective_pos);
        var marker_line = this.lines_in_views[view_idx];
        var end_pos = view.seq2PixScrolled((this.span || 1)  + effective_pos);
        if (end_pos < pix_pos) {
            var pos = end_pos;
            end_pos = pix_pos;
            pix_pos = pos;
        }
        var pix_width = end_pos - pix_pos;
        // small correction for visible offset
        pix_pos += 1;
        end_pos += 1;
        if (!this.span && pix_width < 9) {
            pix_pos = view.seq2PixScrolled(this.seq_pos);
            pix_width = 1;
            marker_line.setStyle("border-right", "none");
        }
        else marker_line.setStyle("border-right", "1px solid " + this.color);
        // Some browsers, WebKit-based in particular, don't like coordinates too far
        // beyond the window. We need to have an ability to hide a marker, that is why we
        // we cull coordinate by +- 10000
        pix_pos = Math.max(-10000, Math.min(10000, pix_pos - SeqView.Marker.getWidth()));
        if (pix_width > 1) {
            end_pos = Math.max(-10000, Math.min(10000, end_pos - SeqView.Marker.getWidth()));
            // Also, we should observe IE limit for width of transparent objects, ca. 4100px
            // Being on the safe side, also compatible with out image size - 4000px
            // Pick the end closest to the mid-window and adjust the other end so that
            // total width is less than 4000.
            if (end_pos - pix_pos > 4000) {
                var mid_view = Math.round(view.getScreenWidth() / 2);
                //console.log("mid_view ", mid_view);
                if (Math.abs(pix_pos - mid_view) < Math.abs(end_pos - mid_view)) {
                    end_pos = Math.max(mid_view + 2000, Math.min(end_pos, pix_pos + 4000));
                    pix_pos = Math.max(pix_pos, end_pos - 4000)
                } else {
                    pix_pos = Math.min(mid_view - 2000, Math.max(pix_pos, end_pos - 4000));
                    end_pos = Math.min(end_pos, pix_pos + 4000);
                }
            }
            pix_width = end_pos - pix_pos;
        }
        marker.setLeft(pix_pos);
        marker_line.setWidth(pix_width);
        // To satisfy IE pre-8
        if (pix_width > 1) {
            marker.setWidth(pix_width + 16);
        }


    },

    //////////////////////////////////////////////////////////////////////////
    // onMouseDown:

    onMouseDown: function(e) {
        //console.log('marker mouse down');
        if (e.button !== 0) { return; }
        if ((this.flags & (SeqView.MarkerFlags.SystemLock | SeqView.MarkerFlags.UserLock)) || Ext.fly(e.getTarget()).hasCls('sv-drag')) {
            this.m_PrevXY = null;
            this.m_CurMarker = null;
            // ignore event for locked marker and pass it to the other handlers
            return;
        }
        //hide tooltip when drugging marker

        if (this.m_ToolTip && !this.m_ToolTip.isPinned()) this.m_ToolTip.hide();
        var marker_node = Ext.fly(e.getTarget()).findParent('div.sv-marker');
        if (!marker_node) { return; }
        if (Ext.isIE) e.stopPropagation(); else e.stopEvent();
        this.initMarkerMove(marker_node.id, e.getXY());
        this.m_DocMouseMove = document.onmousemove;
        this.m_DocMouseUp = document.onmouseup;
        var marker = this;
        document.onmousemove = function(e) { marker.onMouseMove(new Ext.event.Event(e ? e : window.event)); };
        document.onmouseup = function(e) { marker.onMouseUp(new Ext.event.Event(e ? e : window.event)); };
    },

    //////////////////////////////////////////////////////////////////////////
    // onMouseUp:

    onMouseUp: function(e) {
        if (!this.m_PrevXY || !this.m_CurMarker) { return; }
        e.stopEvent();
        this.setCursor(this.m_CurMarker, 'pointer');
        this.m_PrevXY = null;
        this.m_CurMarker = null;
        document.onmousemove = this.m_DocMouseMove;
        document.onmouseup = this.m_DocMouseUp;
        this.m_DocMouseMove = null;
        this.m_DocMouseUp = null;
        this.m_MInfo.updateInfo();

        //updating marker tooltip content
        this.m_ToolTip.updateMarkerToolTipContent(e, this.span, this.seq_pos);
    },

    //////////////////////////////////////////////////////////////////////////
    // onMouseMove:

    onMouseMove: function(e) {
        if (!this.m_PrevXY || !this.m_CurMarker) { return; }
        NCBIGBUtils.ClearBrowserSelection();
        this.moveMarker(e.getXY());
    },

    //////////////////////////////////////////////////////////////////////////
    // onMouseOut:

    onMouseOut: function(e) {
        this.onMouseUp(e);
        //console.log(this);
    },

    //////////////////////////////////////////////////////////////////////////
    // onContextMenu:

    onContextMenu: function(e) {
        e.stopEvent();
    },

    //////////////////////////////////////////////////////////////////////////
    // showSequence:

    showSequence: function() {
        this.m_MInfo.m_App.createTextView(this.seq_pos);
    },

    //////////////////////////////////////////////////////////////////////////
    // centerInView:

    centerInView: function(view) {
        var new_len = view.m_VisLenSeq;
        var new_from = this.seq_pos - view.m_VisLenSeq / 2;
        // Center range of marker, not its beginning if possible
        if (this.span) {
            if (this.span < new_len)
                new_from += Math.floor(this.span / 2);
            else
                new_from += view.m_VisLenSeq / 3;
        }

        if (new_from < 0) new_from = 0;
        if (new_from + new_len > view.m_App.m_SeqLength) new_len = view.m_App.m_SeqLength - new_from;
        if (new_from === view.m_VisFromSeq && new_len === view.m_VisLenSeq) return;
        view.startImageLoading(new_from, new_len);
    },

    //////////////////////////////////////////////////////////////////////////
    // zoomSeqMarker:

    zoomSeqMarker: function(view) {
        view.zoomSeq(this.seq_pos);
    }
};


/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.MarkersInfo
/********************************************************************/


SeqView.MarkersInfo = (function () {

    function constructor(app) {
        this.m_App = app;
        this.m_Dlg = null;
        this.m_AllMarkers = {};
        this.m_AllMarkersRefs = [];
        this.m_MarkerColorIdx = 0;
        this.m_NextMarkerNum = 1;
    };

    constructor.deleteMarker = function(app_idx, marker_num) {
        var app = SeqView.App.findAppByIndex(app_idx);
        if (app && app.m_MarkersInfo) {
             app.m_MarkersInfo.deleteMarker(marker_num);
        }
    };

    constructor.setMarkerPosition = function(app_idx, marker_num) {
        var app = SeqView.App.findAppByIndex(app_idx);
        if (app && app.m_MarkersInfo) {
             app.m_MarkersInfo.setMarkerPosition(marker_num);
        }
    };

    constructor.renameMarker = function(app_idx, marker_num) {
        var app = SeqView.App.findAppByIndex(app_idx);
        if (app && app.m_MarkersInfo) {
             app.m_MarkersInfo.renameMarker(marker_num);
        }
    };

    return constructor;
}) ();

SeqView.MarkersInfo.prototype = {
    getNextColor: function() {
        var sm_MarkerColors = ['green', 'blue', 'red','cyan', '#FF00FF', 'black'];
        var color = sm_MarkerColors[this.m_MarkerColorIdx++];
        if (this.m_MarkerColorIdx > sm_MarkerColors.length-1)
            this.m_MarkerColorIdx = 0;
        return color;
    },

    suggestNextName: function() {
        //var name = this.m_NextMarkerNum;
        var array = this.m_AllMarkersRefs;
        var marker_nums = [];
        for (var i = 0, len = array.length; i < len; i++) {
            var m = array[i];
            if (!m.deleted) {
                var parts = m.marker_name.match(/^Marker (\d+)$/);
                if (parts != null) {
                    marker_nums.push(parts[1]);
                }
            }
        }
        marker_nums.sort(function(a,b){return a-b});
        var next_num = 1;
        for (var i = 0, len = marker_nums.length; i < len; i++) {
            if (next_num != marker_nums[i]) break;
            next_num++;
        }
        return "Marker " + next_num;
 },

//////////////////////////////////////////////////////////////////////////
// showDlg:

    showDlg: function( view ) {
        var active_marker_count = 0;
        this.forEachMarker(function() { active_marker_count++; });

        if (active_marker_count == 0) {
            this.MarkerDlg(view);
            return;
        }
        //if (!this.m_Dlg) {
            this.m_Dlg = new Ext.Window({
                layout:'fit',
                title:'Markers',
                app: view.m_App,
                minWidth:480, width:600, height:350,
                cls: 'SeqViewerApp',
                constrain: true, 
                collapsible: true,
                listeners: {
                    close: {scope: this, fn: function() {
                        this.m_App.un('origin_changed', this.onOriginChanged, this);
                        this.m_App.un('strand_changed', this.onStrandChanged, this);
                    }}
                },
                tbar:[
                    { text: 'Add Marker', iconCls:'xsv-marker-add', tooltip:'Add New Marker', scope:this,
                        handler:function() { this.MarkerDlg(view); }
                    },
                    '->',
                    { text:'Remove all markers', iconCls:'xsv-clear_markers', scope:this,
                        handler:function() {
                            Ext.MessageBox.confirm( 'Confirm', 'Clear all markers?', function(btn) {
                                if (btn!='yes') return;
                                this.reset();
                            },
                            this);
                        }
                    }
                ],
                closeAction: 'destroy',
                items:[{
                    xtype:'panel', autoScroll: true, itemId: 'infoPanel', layout:'fit'
                }],
                buttons: [
                    {text: 'Download', scope: this, handler: function(button) {
                        this.m_App.pingClick('7-2-1');
                        var form = Ext.DomHelper.append(document.body, 
                            { tag : 'form', method : 'post',  action : this.getMarkerInfoURL(true) });
                        document.body.appendChild(form);
                        form.submit();
                        document.body.removeChild(form);
                    }},

                    {text: 'Close', scope: this, handler: function() { this.m_Dlg.close(); this.m_Dlg = null;} }
                ]
            });
        //}
        this.m_App.on('origin_changed', this.onOriginChanged, this);
        this.m_App.on('strand_changed', this.onStrandChanged, this);

    var app = this.m_App;
    this.updateInfo();    
    this.m_Dlg.on( 'close', function() { app.resizeIFrame(); app.m_DialogShown = false; } );

    app.resizeIFrame( 400 );
    app.m_DialogShown = true;
    this.m_Dlg.show();
    },

    onOriginChanged: function(app) { this.updateInfo(); },
    onStrandChanged: function(app) { if (app.getOrigin() != 0) this.updateInfo(); },

    addMarker: function(range, name, lock, color, flags) {
        var newMarker = new SeqView.Marker(this, range, name, lock, color, flags);
        this.m_AllMarkersRefs.push(newMarker);
        this.updateInfo();
        return newMarker;
    },

//////////////////////////////////////////////////////////////////////////
// addMarkerByData:

    addMarkerByData: function(data) {
        this.addMarker(data[0], data[1], data[2], data[3], data[4]);
    },

    MarkerDlg: function(view, pos, oldMarker){
        var setMarkerDlg;

        if (view) { // pos is X coordinate on the image
            pos = pos || Math.floor(view.getScreenWidth()/2);
            pos = this.m_App.posToLocalDisplay(view.pix2Seq(pos - view.m_ScrollPix));
        } else { // if pos set it's position on the sequence
            pos = pos || this.m_App.posToLocalDisplay(Math.floor(this.m_App.m_SeqLength / 2));
            view = this;
        }

        setMarkerDlg = new Ext.Window({
            modal:true, title:'Set New Marker',
            app: view.m_App,
            width:270,
            collapsible:false, resizable:false, closeAction:'close', iconCls:'xsv-markers',
            itemId: 'setMarkerDlg',
            border: false,
            buttonAlign:'center',
            viewModel: { data: {
                color: oldMarker ? oldMarker.color : this.getNextColor(),
                lock: false,
                name: oldMarker ? oldMarker.visible_name : this.suggestNextName(),
                position: pos
            }},
            items:[
               {xtype: 'form', frame: true,
                items:[
                   {xtype:'textfield', width:'100%', fieldLabel: 'Name', hidden: oldMarker != null, bind: '{name}'},
                   {xtype:'textfield', width:'100%', fieldLabel: 'Pos/Range',  bind: '{position}'},
                   {xtype: 'container', layout: 'hbox',  height: 22,
                    items: [
                        {xtype:'displayfield', fieldLabel: 'Color'},
                        {xtype:'component', width: 22, height: 22, bind: {style: {backgroundColor:  '{color}'}}},
                        {xtype:'button',
                            menu: { xtype: 'colormenu',
                                handler: function(picker, choice) {
                                    picker.up('#setMarkerDlg').getViewModel().setData({color:'#' + choice});
                                } 
                        }}]
                   },
                   {xtype:'checkbox', height:22, labelSeparator:'', boxLabel:'Lock Marker', hidden: oldMarker != null, bind: '{lock}'}
                ]}
            ],
            buttons:[
               {text:'OK', scope:this, handler: function() {
                  var data = setMarkerDlg.getViewModel().getData();
                  if (!data.name ||  data.name.length > 50) {
                      Ext.MessageBox.show({ title: 'Set New Marker', msg: 'Invalid marker name.',
                                            buttons: Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                      return;
                  }
                  this.m_App.handlePos(data.position.toString(), {
                      ask_user: true,
                      success: function(pos_range, options) {
                          setMarkerDlg.close();
                          if (oldMarker) this.deleteMarker(oldMarker.marker_num);
                          this.addMarker(pos_range, data.name, data.lock, data.color);
                      },
                      failure: function(err_msg, options) {
                          Ext.MessageBox.show({title: 'Set New Marker', msg: err_msg,
                              buttons: Ext.MessageBox.OK, icon: Ext.MessageBox.ERROR });
                      },
                      scope: this
                  });
               }},
               {text:'Cancel', handler: function() {setMarkerDlg.close(); } }
             ]
        });
        setMarkerDlg.show();
    },


    getMarkerInfoURL: function(for_download, format) {
        var url = this.m_App.m_CGIs.ObjCoords + '?objcoords=1&id=' + this.m_App.GI 
        var positions = [];
        var names = [];
        this.forEachMarker(function(m) {
            positions.push(m.seq_pos);
            names.push(m.visible_name);
            if (m.span) {
                positions.push(m.seq_pos + m.span - 1);
                names.push(m.visible_name);
            }
        });
        if (for_download) {
            url += '&download=true&fmt=' + (format ? format : 'csv');
            for (var i = 0; i < positions.length; i++) {
                url += '&pos=' + positions[i] + '&name=' + names[i];
            }
        } else {
            var sorted = positions.sort(function(a,b) { return a - b;});
            for (var i = 0, prev = sorted[0] - 1; i < sorted.length; i++) {
                if (sorted[i] != prev) url += '&pos=' + sorted[i];
                prev = sorted[i];
            }
        }
        return url;
    },


    parseObjCoords: function(obj_coords_json, marker) {
        var data = obj_coords_json;
        var pos  = marker.seq_pos;
        var end_pos;
        if (marker.span)
            end_pos = pos + marker.span - 1;
        var html = '';
        var useHGVS = this.m_App.m_ViewParams.organism == "Homo sapiens";
        var cur_pos = -1;
        for (var i = 0; i != data.length; i++) {
            var row_data = data[i];
            var marker_pos = row_data['marker_pos'];
            if (marker_pos != pos  &&  marker_pos != end_pos) continue;
            if (html.length == 0) {
                html += '<table class="xsv-markerinfo" width="98%">'
                html += '<tr><th><b>Name:</b> ' + marker.visible_name;
                if (!marker.nvMarker) html += ' (<a href="#" onClick="SeqView.MarkersInfo.renameMarker('+this.m_App.m_Idx+','+marker.marker_num+');">Edit</a>)</th>';
                html += '<th width="75%" colspan=' + (useHGVS ? 4 : 3) + '>Position: ' + marker.displayPos();
                if (!marker.nvMarker) html += ' (<a href="#" onClick="SeqView.MarkersInfo.setMarkerPosition('+this.m_App.m_Idx+','+marker.marker_num+');">Edit</a>)';
                html += '<a href="#" onClick="SeqView.MarkersInfo.deleteMarker('+this.m_App.m_Idx+','+marker.marker_num+');" class="xsv-marker_remove_link">Remove</a>';
                html += '</th></tr>';
                if (useHGVS)
                    html +='<tr><th width="25%">Accession/Locus tag</th><th width="10%">Location</th><th width="10%">Relative to</th><th width="30%"><a href="http://www.hgvs.org/mutnomen" target="_blank">HGVS Name</a></th><th width="25%">Sequence</th></tr>';
                else
                    html +='<tr><th width="25%">Accession/Locus tag</th><th width="20%">Location</th><th width="10%">Relative to</th><th width="45%">Sequence</th></tr>';
            }
            html += this.addMarkersRow(row_data, useHGVS, false, cur_pos != marker_pos);
            if (cur_pos != marker_pos  &&  this.m_App.m_Origin)
                html += this.addMarkersRow(row_data, useHGVS, true);
            cur_pos = marker_pos;
        }
        html += '</table>';
        return html;
    },

    addMarkersRow: function(data, useHGVS, use_local, use_separator) {
        var pos_mapped = '';
        if (typeof(data['pos_mapped']) != 'undefined') {
            pos_mapped = data['pos_mapped'];
            if (use_local) {
                pos_mapped = this.m_App.posToLocal(pos_mapped);
            } else {
                if (pos_mapped >= 0) pos_mapped += 1;
            }
        }
        if (!useHGVS && pos_mapped.length == 0) return '';

        var html = use_separator ? '<tr class="sv-rowgroup">' : '<tr>';
        html += '<td>' + data['title'] + '</td>'; // Accession
        html += '<td>' + pos_mapped + '</td>'; // Location
        var hgvs = data['hgvs_position'];
        var rel_to = "Seq start";
        if (use_local) {
            rel_to = "Current origin";
        }
        else if (hgvs  &&  hgvs.charAt(0) == 'c') {
            rel_to = "CDS start";
        }
        html += '<td>' + rel_to + '</td>'; // Relative to
        if (useHGVS) {
            if (hgvs  &&  data['title']) {
                hgvs = data['title'] + ':' + hgvs;
            }
            html += '<td>' + hgvs + '</td>';
        }
        html += '<td><span style="font-family:monospace">' + data['sequence']+'</span></td>';
        html += '</tr>';
        return html;
    },


    updateInfo: function() {
        if (!this.m_Dlg) return;
        
        var loadCallback = function(data, txt, rq){
            var formated_html = '';
            this.forEachMarker(function(m) { formated_html += this.parseObjCoords(data, m); }, this);
            if (formated_html.length == 0) { // no markers to display
                formated_html = '<div style="color:gray;margin-top:10px;" align="center">No Markers To Display</div>';
            }
            this.m_Dlg.down('#infoPanel').update(formated_html);
        }
        SeqView.App.simpleAjaxRequest({url: this.getMarkerInfoURL() + '&appname=' + this.m_App.m_AppName, context: this,
        success: loadCallback, error: loadCallback});
    },

    deleteMarker: function(num) {
        var i = this.forEachMarker(function(m) {
            if (m.marker_num == num) {
                m.deleteMarker();
                return false;
            }
        });
        if (i)
            this.updateInfo();
    },

    setMarkerPosition: function(mNum) {
        var marker = null;
        this.forEachMarker(function(tmp) { if (tmp.marker_num == mNum) marker = tmp; });
        if (marker) this.MarkerDlg(null, marker.displayLocalPos(), marker);
    },

//////////////////////////////////////////////////////////////////////////
// renameMarker:

    renameMarker: function(num) {
        var m = null;
        this.forEachMarker(function(tmp) {
            if (tmp.marker_num == num) {
                m = tmp;
            }
        });
        if (!m) return; // marker not found

        Ext.MessageBox.prompt('Marker', 'Please enter new marker name:', function(btn, text) {
            if (btn!='ok' || text.length==0) return;
            m.rename(text);
        }, this, false, m.marker_name );
    },

//////////////////////////////////////////////////////////////////////////
// reset:

    reset: function() {
        for(var m in this.m_AllMarkers) {
            if (this.m_AllMarkers[m]) {
                this.m_AllMarkers[m].deleteMarker();
            } // delete markers
        }
        this.updateInfo();
        this.m_AllMarkersRefs = [];
        this.m_MarkerColorIdx = 0;
        this.m_NextMarkerNum = 1;
    },

//////////////////////////////////////////////////////////////////////////
// updateMarkersSize:

    updateMarkersSize: function(view) {
        this.forEachMarker(function(m) { m.updateMarkerSize(view); });
    },

//////////////////////////////////////////////////////////////////////////
// updateMarkersPos:

    updateMarkersPos: function(view) {
        this.forEachMarker(function(m) {
            m.updateMarkerPos(view);
            //console.log('update MarkerPosition');
            //update marker tooltip range text content
            m.m_ToolTip.updateMarkerToolTipContent(null,this.span,this.seq_pos);
        });
    },

//////////////////////////////////////////////////////////////////////////
// scrollMarkers:

    scrollMarkers: function(view, delta) {
        this.forEachMarker(function(m) { m.scrollPix(view, delta); });
    },

//////////////////////////////////////////////////////////////////////////
// forEachMarker: call fn in scope for each marker, if return of fn is false
//   return number of marker plus one for this function and cancel execution,
//   otherwise return 0.
    forEachMarker: function(fn, scope) {
        var array = this.m_AllMarkersRefs;
        for(var i = 0, len = array.length; i < len; i++){
            var m = array[i];
            if (!m.deleted) {
                if (fn.call(scope || m, m, i, array) === false) { return (i+1); }
            }
        }
        return 0;
    },

//////////////////////////////////////////////////////////////////////////
// findMarker:

    findMarker: function(m_name) {
        return this.m_AllMarkers[m_name];
    },

//////////////////////////////////////////////////////////////////////////
// findMarkerByName:

    findMarkerByName: function(name) {
        var array = this.m_AllMarkersRefs;
        for(var i = 0, len = array.length; i < len; i++){
            var m = array[i];
            if (!m.deleted && m.marker_name == name) {
                return m;
            }
        }
        return null;
    },

//////////////////////////////////////////////////////////////////////////
// findMarkerByPos:

    findMarkerByPos: function( pos ){
        var array = this.m_AllMarkersRefs;
        for(var i = 0, len = array.length; i < len; i++){
            var m = array[i];
            if( !m.deleted && m.seq_pos == pos ){
                return m;
            }
        }
        return null;
    },

//////////////////////////////////////////////////////////////////////////
// hasMarkers:

    hasMarkers: function() {
        return this.m_AllMarkersRefs.length > 0;
    },

    parseMarkersURL: function(markers, val) {

    if( !val || val.length == 0 ){
        return;
    }

        var parts = cbSplit(val, ',');
        while (parts.length) {
            var part = parts.splice(0,1)[0];
            while (parts.length && part.charAt(part.length-1) == '\\') {
                part = part.substr(0, part.length-1) + ',' + parts.splice(0,1)[0];
            }

            var lock = false;
            var lastChar = part.charAt(part.length-1);
            if (lastChar === '!') {
                lock = true;
                part = part.substr(0, part.length-1);
            }

            var info   = part.split('|');
            var coords = info.splice(0,1)[0].split(':');
            var name   = info.splice(0,1)[0];
            while (info.length && name.charAt(name.length-1) == '\\') {
                name = name.substr(0, name.length-1) + '|' + info.splice(0,1)[0];
            }
            name = NCBIGBUtils.unescapeName(name);
            var color  = info.splice(0,1)[0];
            var flags  = info.length > 0 ? parseInt(info.splice(0,1)[0]) : 0;
            if (lock) flags = flags | SeqView.MarkerFlags.UserLock;

            if( isNaN(coords[0]) ) continue;

            var range;
            var marker_pos = NCBIGBUtils.stringToNum(coords[0])-1;
            if( coords.length > 1 ){
                if( isNaN(coords[1]) ) continue;

                range = [marker_pos, NCBIGBUtils.stringToNum(coords[1])-1]
            } else {
                range = [marker_pos];
            }

            var colors = { red: "ff0000", green: "008000", blue: "0000ff",
                           black: "000000", cyan: "00ffff" };
            if (color && color[0] != '#') {
                if (/^[0-9a-fA-F]+$/.test(color)) {
                    color = '#' + color;
                } else if( colors[color] ){
                    color = '#' + colors[color];
                } else {
                    color = '#' + colors['green'];
                }
            }
            markers.push([range, name, flags, color]);
        }
    },

    parseMarkersOldStyle: function(markers, marker_ranges, marker_names) {
        for (var m = 0; m < marker_ranges.length; m++) {
            var val = marker_ranges[m];
            var flags = 0;
            if (val.indexOf('!') != -1) {
                flags = flags | SeqView.MarkerFlags.UserLock;
                val = val.replace(/!/, '');
            }
            var parts = val.split(':');
            var marker_pos = NCBIGBUtils.stringToNum(parts[0])-1;
            var range;
            if (parts.length > 1) {
                range = [marker_pos, NCBIGBUtils.stringToNum(parts[1])-1]
            } else {
                range = [marker_pos];
            }
            markers.push([range, marker_names[m], flags]);
        }
    },

//////////////////////////////////////////////////////////////////////////
// getMarkersData:
// if less_strict == true -- only symbols '|' and ',' will be "escaped"
    getMarkersData: function(less_strict, forPDF) {
        var markers = "";

        this.forEachMarker(function(m) {
            if (!m.seq_pos || isNaN(m.seq_pos) || m.nvMarker) return;
            if (markers.length > 0) markers += ',';

            markers += m.displayPos(forPDF);
            markers += "|" + NCBIGBUtils.escapeName(m.marker_name, less_strict ? /([\|,])/ : false);
            markers += "|" + ((m.color.charAt(0) === '#') ? m.color.slice(1) : m.color);
            if (m.flags) markers += "|" + m.flags;
        });

        return markers;
    },

//////////////////////////////////////////////////////////////////////////
// getMarkersURL:

    getMarkersURL: function() {

        var mk_data = this.getMarkersData();
        if( mk_data && mk_data.length > 0 ){
            mk_data = '&mk=' + mk_data;
        }

        var dlg_url = '';
        if( this.m_Dlg && this.m_Dlg.isVisible() ){
            dlg_url = '&vm=true';
        }

        return mk_data + dlg_url;
    }
};

Ext.define('SeqView.MarkerToolTip', {
    extend: 'NCBIGBObject.ToolTip',
    width: 180,
    anchor: 'top',
    initComponent: function() {
        this.title = '&nbsp;&nbsp;&nbsp;' + (this.selection.nvMarker ? 'Navigation Marker' : this.selection.marker_name);
        this.textId = Ext.id();
        this.setMarkerRange(this.selection.seq_pos, this.selection.span);
        this.callParent(arguments);

        var items = [{xtype: 'tbtext', text: this.rangeText, id: this.textId, textStyle: 'font-weight:bold;'},
                    {text: (this.selection.nvMarker ? 'Set New Marker at Position' : 'Rename...'), iconCls:'xsv-markers',
                    handler: function() {
                        if (!this.selection.nvMarker) { 
                            this.selection.m_MInfo.renameMarker(this.selection.marker_num);
                            this.pingClick('7-N');

                        } else {
                            this.selection.m_MInfo.MarkerDlg(null, this.selection.seq_pos + 1);                            
                            this.pingClick('7-R');
                        }
                    }}];
        if (!this.selection.nvMarker) {
            items.push(
                {text: 'Modify Position/Range/Color', iconCls: 'xsv-zoom_range',
                handler: function() {
                    this.pingClick('7-0');
                    if (!this.isPinned()) this.hide();
                    if (!(this.selection.flags & (SeqView.MarkerFlags.SystemLock | SeqView.MarkerFlags.UserLock))) {
                        this.selection.m_MInfo.setMarkerPosition(this.selection.marker_num);
                    } else 
                        Ext.Msg.alert('Information',
                            (this.selection.flags & SeqView.MarkerFlags.UserLock)
                            ? 'To modify range please first unlock the marker.'
                            : "The object is locked, can't modify range");
            }});
        }
        items.push(
            {text: 'Zoom To Sequence At Marker', iconCls:'xsv-zoom_seq',
            handler: function() {
                this.pingClick('7-1');
                this.selection.zoomSeqMarker(this.gview);
                if (!this.isPinned()) this.hide();
            }},
            {text: 'Marker Details', iconCls: 'xsv-markers',
            handler: function() {
                this.pingClick('7-2');
                if (!this.isPinned()) this.hide();
                this.selection.m_MInfo.showDlg(this.gview);
        }});
       

        var tooltip = this;
        var selection = this.selection;
        var app = this.selection.m_MInfo.m_App;
        if (!app.m_Embedded) {
            var vNum = 0;
            app.forEachView(function(view) {
                if (view.isGraphic() && !(vNum++)) {
                    items.push({
                         text: 'Reveal in Graphic View', iconCls: 'goto_marker', scope: this,
                         handler: function() {
                             this.centerInView(view);
                             this.m_ToolTip.pingClick('7-Graph');
                         }
                    });
                }
            }, selection );
        }
        items.push(
             {text: 'Reveal in Sequence View', iconCls: 'xsv-new_fasta', scope: selection,
             handler: function() {
                 this.showSequence();
                 this.m_ToolTip.pingClick('7-Seq');
             }
        });
        if (!(this.selection.flags & SeqView.MarkerFlags.SystemLock) && !this.selection.nvMarker) {
            items.push(
                {text: 'Lock/Unlock Marker', iconCls: 'xsv-marker_locked',
                handler: function() {
                    this.selection.lockMarker();
                    this.pingClick('7-3');
                }
            });
        }
        items.push(
            {text: 'Remove Marker', iconCls: 'xsv-marker-remove',
            handler: function() {
                if (!this.selection.m_ToolTip.pinned) {
                    this.selection.deleteMarker();
                    this.selection.m_MInfo.updateInfo();
                    this.pingClick('7-4');
                } else {
                    Ext.Msg.alert('Status', 'Can\'t remove marker with pinned tooltip.');
                }
            }
        });

        var scope = this;
        Ext.each(items, function(item, ix){
            item.scope = item.scope || scope; 
            item.style = {margin: '0px', padding: '0px'};
            item.x = 0;
            item.y = 20 * ix;
        });        
        this.toolbar = new Ext.Toolbar({layout: 'absolute', target: this.target, items: items, border: false,
            height: items.length * 20});
       
        this.add(this.toolbar);
        this.addTool([
            { type: 'search', hidden: true, scope: this,
                callback: function() {
                    if (!this.isPinned()) this.hide();
                    this.centerMarker();
                }           
            },{ type: 'magnify', scope: this,
                callback: function() {
                    this.selection.zoomSeqMarker(this.gview);
                    if (!this.isPinned()) this.hide();
                }
            }]);
    },

    centerMarker: function () {
        var seq_end = this.selection.seq_pos + this.selection.span;
        var range_l = this.selection.seq_pos;
        var range_r = seq_end;
        var pos_l = this.gview.seq2Pix(range_l);
        var pos_r = this.gview.seq2Pix(range_r);
        var view_start = -this.gview.m_ScrollPix;
        var valid = true;
        if (pos_l < 0 || pos_r > 4000) {
            valid = false;
        }
        var jj = this.gview.m_FromCgi.areas;
        if (valid) {
            if (pos_r - 10 < view_start) {
                this.gview.scrollViewTo(-(pos_l-10), SeqView.PAN_RIGHT);
            } else  if (pos_l+10 > (this.gview.getScreenWidth() + view_start)) {
                var new_pos = pos_r - this.gview.getScreenWidth();
                this.gview.scrollViewTo(-(new_pos+10), SeqView.PAN_LEFT);
            }
        } else {
           this.selection.zoomSeqMarker(this.gview);
        }
    },
    
    setMarkerRange: function(seq_pos, span) {
        this.adjusted_seq_pos = this.gview.m_App.posToLocal(seq_pos);
        this.adjusted_seq_end = null;
        if (span) {
            this.adjusted_seq_end = this.gview.m_App.posToLocal(seq_pos + span - 1);
            this.rangeText = '<b>Range: ';
            if (this.gview.m_App.getFlip()) {
                this.rangeText += this.adjusted_seq_end + ' .. ' + this.adjusted_seq_pos;
            } else {
                this.rangeText += this.adjusted_seq_pos + ' .. ' + this.adjusted_seq_end;
            }
            this.rangeText += '</b>';
        } else
            this.rangeText = '<b>Position: ' + this.adjusted_seq_pos + '</b>';
    },
    
    updateMarkerToolTipContent: function(e,span,seq_pos) {
        this.setMarkerRange(seq_pos, span);
           var el = Ext.getCmp(this.textId);
           el.setText(this.rangeText);
    },

    isPinned: function() {
        return this.pinned;
    },

    onTargetOut: function(e){
        if(this.el && !this.disabled && !e.within(this.target.dom, true)){
            this.el.applyStyles('outline-style: none');
        }
        this.callParent(arguments);
    },
    onTargetOver: function(e){
        if (this.isPinned()) {
            if (!this.isVisible()) {
                this.mouseOffset = [0, 0];
                this.show();
            }
            this.el.applyStyles('outline-style: solid;outline-color: red;outline-width:thin');
        }
        else this.callParent(arguments);
//        this.ping({"jsevent":"mouseover", "sv-event":"Marker_ToolTip"});
    },
    onShow: function() {
        if( this.gview.m_ContextMenu && this.gview.m_ContextMenu.isVisible() ){
                this.gview.m_ContextMenu.hide();
        }
        this.callParent(arguments);
    },
    ping: function(a) { this.selection.m_MInfo.m_App.ping(a); },
    pingClick: function(a, e) { this.selection.m_MInfo.m_App.pingClick(a, e); }

});
/*  $Id: view.js 37983 2017-03-09 23:18:17Z borodine $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vlad Lebedev, Maxim Didenko
 *
 * File Description:
 *
 */

SeqView.colorBarTPL = new Ext.Template('<div id="view_color_{idx}" style="width:15px;height:15px;background-color:#{color};border:1px black solid;">&nbsp;</div>');

Ext.define('SeqView.View', {
    constructor: function(type, app) {
        this.clear();
        this.m_Type = type;
        this.m_App = app;
        this.m_ReqNum = 0;
        SeqView.m_NextViewIdx = SeqView.m_NextViewIdx || 0;
        this.m_Idx = SeqView.m_NextViewIdx++;
        this.m_View = null;
    },

    getSpacerHeight: function() { return 4; },
    getWidth: function() { return this.m_Width; },
    getHeight: function() { return this.m_Height; },
    isPanorama: function() { return false; },
    isGraphic: function() { return false; },
    isAlignment: function() { return false; },
    getMostRightPix: function() { return 0; },
    getScreenWidth: function() { return this.m_App.getWidth() - (this.isPanorama() ? 2 : 0); },
    getXY: function() { return this.m_View.getEl().getXY(); },
    ping: function(a) { this.m_App.ping(a); },
    pingClick: function(a, e) { this.m_App.pingClick(a, e); },
    updateTitle: function() {},
    refresh: function() {},
    isLoading: function() { return this.m_Loading; },
    moveTo: function(vis_from, vis_len) {},

    createChooseColorBtn: function() {
        var m_ncbi_app = this.m_App.m_ncbi_app;
        if (!this.m_Color || this.m_Color.length == 0) {
            var cp = new Ext.ColorPalette();
            this.m_Color = cp.colors[Math.floor( Math.random() * cp.colors.length)];
        }
        return {
                text:SeqView.colorBarTPL.apply({idx:this.m_Idx, color: this.m_Color}),
                tooltip:'View Color',
                handler: function() {this.pingClick('1-2-0');},
                menu:new Ext.menu.ColorPicker({listeners: {'select': function(cm, color) {
                    Ext.get('view_color_'+this.m_Idx).setStyle('background-color', '#'+color);
                    this.m_Color = color;
                    if (this.m_Locator) { this.m_Locator.setColor(color); }
                    this.m_App.reCreateReflections();
                },scope:this}})
        };
    },

    createCloseTBar: function() {
        return  {
            id:'close', qtip:'Close View', scope:this,
            handler:function(e, target, panel) {
                this.m_App.viewIsClosing(this);
                var view = panel.view;
                if(view) { view.remove(); }
            }
        };
    },

    clear: function() {
        this.m_Loading =  false;
        this.m_TopOffset = this.m_ScrollPix = this.m_Width = this.m_Height = 0;
        this.m_Theme = this.m_DivId = this.m_FromCgi = this.m_View = this.m_Color = this.m_Spacer = this.m_Locator = null;
    },

    destroy: function() {
        if (this.m_Locator) {
            this.m_Locator.remove();
        }
        if(this.m_View) {
            if (this.m_Spacer) {
                this.m_View.ownerCt.remove(this.m_Spacer, true);
            }
            this.m_View.ownerCt.remove(this.m_View, true);
        }
        this.clear();
    },

    remove: function() {
        this.destroy();
        this.m_App.removeView(this);
    },
    setGeneMode: function(){},

    addURLParams: function(params) {
        var app = this.m_App;
        params = params || {};
        if (app.m_NAA) params.naa = app.m_NAA;
        if (app.m_BamPath) params.bam_path = app.m_BamPath;
        if (app.m_SRZ) params.srz = app.m_SRZ;
        if (app.m_Key) params.key = app.m_Key;
        if (app.m_Origin != 0 ) params.origin = app.m_Origin;
        if (app.m_DepthLimit) params.depthlimit = app.m_DepthLimit;
        if (app.m_AssmContext) params.assm_context = app.m_AssmContext;
        if (app.m_highlightsColor) params.highlights_color = app.m_highlightsColor;

        if (app.m_Config) {
            Ext.apply(params, app.m_Config.visualOptionsUrl())
        }
        if (this.isGraphic()) {
            var sm = this.m_slimMode ? ',show_title:false' : '';
            var noParallel =  (params.select && params.select.indexOf(';gi') > 0) && (app.m_hiddenOptions.indexOf('always_parallel') < 0);
            if (params.forPDF || noParallel /*(params.select && params.select.indexOf(';gi') > 0)*/)
                params.tracks = [SeqView.TM.tracksArrayToString(app.getActiveTracks(), true, params.forPDF, sm)];
            else
                params.tracks = SeqView.TM.tracksToArrayOfStrings(app.getActiveTracks(), sm);
            if (!params.tracks.length || !params.tracks[0]) params.tracks[0] = '[key:ruler]';
            if (app.m_ViewLabels) params.markers = app.m_ViewLabels;
            if (app.m_SnpFilter) params.snp_filter = app.m_SnpFilter;
        }
        return params;
    }    
});
/*  $Id: panorama.js 38272 2017-04-19 20:24:16Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko
 * File Description:
 */
 

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Panorama 
/********************************************************************/


Ext.define('SeqView.Panorama', {
    extend: 'SeqView.View',
    m_PrevXY: null,
    m_CurLocator: null,
    m_ActionNeeded: false,

    constructor: function(app) {
        //SeqView.Panorama.superclass.constructor.apply(this, ['panorama', app]);
        this.callParent(['panorama', app]);
        this.m_TopOffset = 18; //sm_HolderSize
        var pID = Ext.id();
        this.m_DivId = 'panorama_id_' + pID;
        this.m_selectionDivId = 'panorama_selection_' + pID;
        this.m_View = this.m_App.addView({ 
            html:'<div id="'+this.m_DivId+'" class="panorama_div"><div id="'+ this.m_selectionDivId + '" class="panoramaSelection" ></div><div id="pan-holder" class="pan-holder"/>'
        });
        Ext.get(this.m_DivId).on({ scope: this,
            mousedown:    this.onMouseDown,
            mouseup:      this.onMouseUp,
            mousemove:    this.onMouseMove,
            touchstart:   this.onMouseDown,
            touchend:     this.onMouseUp,
            touchmove:    this.onMouseMove,
            contextmenu:  this.onContextMenu,
            mouseout:   function() {this.m_vLine.setStyle('display', 'none');}
        });
        //if  embedded=panorama then do not allow selection
        if (!this.m_App.m_AllViewParams.match("embedded=panorama")) 
            this.m_PanoramaSelection = new SeqView.PanoramaSelection(this.m_DivId,this.m_selectionDivId, [0,0]);
        
    },

    isPanorama: function() { return true; },

    loadImage: function() {
        this.m_Width = this.getScreenWidth();
        if (this.m_Width <= 0) {
            // the view is hidden, don't do anything now.
            // the client will have to call refresh method to load the image
            return;
        }
        this.m_Loading = true; // start loading
      
        var url = this.m_App.m_CGIs.Panorama;
        var params = {theme:'NCBI Overview', id: this.m_App.GI, width: this.m_Width};
        Ext.apply(params, this.m_App.m_GraphicExtraParams);
        this.addURLParams(params);
      
        this.m_App.AjaxRequest({url: url, context: this, data: params,
                    success:this.checkJobStatus, error:this.loadFailure});      
    }, 
    
    checkJobStatus: function(data, text, res) {
        var from_cgi = SeqView.decode(data);
        if (from_cgi.job_status) {
            if (from_cgi.job_status == 'failed') {
                this.loadFailure(null, from_cgi.error_message);
            } else if(from_cgi.job_status == 'canceled') {
                this.loadFailure(null, 'Job canceled');
            } else {
                var url = this.m_App.m_CGIs.Panorama + '?job_key=' + from_cgi.job_id
                Ext.defer(this.m_App.AjaxRequest, 2000,this,[{url: url, context: this,
                        success: this.checkJobStatus, error: this.loadFailure}]);
            }
        } else {
            if (from_cgi.success === false) {
                this.loadFailure(null, from_cgi.msg);
            } else {
                // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
                // This is a way to provide reliable URL resolution for embedding. SV-1760
                if (from_cgi.img_url && from_cgi.img_url.charAt(0) == '?') {
                    from_cgi.img_url = this.m_App.m_CGIs.NetCache + from_cgi.img_url;
                }
                this.m_FromCgi = from_cgi;
                this.m_Height = this.m_FromCgi.img_height + this.m_TopOffset + 1;
    
                var the_div = Ext.get(this.m_DivId);
                var img_el = the_div.first('img');
                if (img_el) {
                    var d = Ext.getDom(img_el);
                    d.src = from_cgi.img_url;
                } else {
                    var tpl = new Ext.Template('<img src="{img_url}">'),
                    img_el = tpl.append(the_div, from_cgi,true);
                }
                the_div.setHeight(this.m_Height); 
//                this.m_View.setHeight(this.m_Height + 2);
                this.m_vLine = the_div.appendChild(new Ext.Element(document.createElement('div')));
                this.m_vLine.set({style: 'position:absolute;border:1px solid lightgray;border-right-width:0; top:0; width:0px; height:'
                    + this.m_Height + 'px; display: none;'});
                         
                this.m_App.updateMarkersSize(this); // update markers
                this.m_App.updateMarkersPos(this); // update markers
                this.m_Loading = false; // loaded
                this.m_App.notifyViewLoaded(this);

                if( !this.m_App.m_DialogShown ){
                    this.m_App.resizeIFrame();
                }
            }
        }
    },
    loadFailure: function(data, text, res) {
        Ext.MessageBox.show({title:'Image loading error', msg: text, buttons: Ext.MessageBox.OK, icon: Ext.MessageBox.INFO});
    },


    toPix: function(x) { // sequence 2 screen (panorama)
        return x * this.m_Width / this.m_App.m_SeqLength;
    },

    toSeq: function(x) { 
        return Math.round(this.m_App.m_SeqLength * x / this.m_Width);
    },
    
    seq2Pix: function(seq_pos) {
        return this.toPix(seq_pos);
    },
    
    seq2PixScrolled: function(seq_pos) {
        // For panorama m_ScrollPix is always zero, but anyway
        return this.seq2Pix(seq_pos) + this.m_ScrollPix;
    },

    pix2Seq: function(pix_pos) { // screen 2 sequence (gview)
        return this.toSeq(pix_pos);
    },

    refresh: function() {
        this.loadImage();
    },

    getMostRightPix: function() {
        return this.getWidth();
    },
    
    onMouseDown: function(e) {
        if (this.m_ContextMenu) this.m_ContextMenu.destroy();
        if (e.type == 'mousedown' && e.button) return
        var tID = e.target.id;
        if (this.m_PanoramaSelection && tID.search('scroller') < 0 && tID.search('resizer') < 0){
            this.m_ResizeAction = true;
            this.m_XY = e.getXY();
            if (e.type != 'mousedown') {
                this.m_deferredContext = Ext.defer(this.onContextMenu, 2000, this, null);
            } else this.m_deferredContext = 0;
            this.m_XFinal_Selection = this.m_XY[0] + 1;
            this.m_PanoramaSelection.resize(this.m_XY[0], this.m_XFinal_Selection);
            var el =  Ext.get(this.m_selectionDivId);
           if (el) el.applyStyles("display:inline;");
           e.stopEvent();
        }
    },

    onMouseMove: function(e) {
        var magicSensivity = 5; 
        var x = e.getX();
        if (!this.m_ResizeAction || !this.m_PanoramaSelection) {
            this.m_vLine.setStyle('display', 'block');
            this.m_vLine.setX(x);
            return;
        }
        this.m_vLine.setStyle('display', 'none');

        if (this.m_deferredContext) {
            if (Math.abs(x - this.m_XY[0]) <= magicSensivity) return;   
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        this.m_XFinal_Selection = x;
        this.m_PanoramaSelection.resize(this.m_XY[0], x);
        this.m_Selection = true;
        e.stopEvent();
    },

    onMouseUp: function(e) {
        var el =  Ext.get(this.m_selectionDivId);
        if (!el || !el.isVisible() || !this.m_PanoramaSelection) return;
        if (this.m_deferredContext) {   
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        this.m_ResizeAction = false;
        if (Math.abs(this.m_XFinal_Selection - this.m_XY[0]) < 5) {
            el.applyStyles("display:none;");
            return;
        } else {
            if(this.m_App.m_Views[1]){
                var m_Locator = this.m_App.m_Views[1].m_Locator;
                m_Locator.setLeft(el.getLeft() - Ext.get(this.m_DivId).getLeft());
                m_Locator.setWidth(el.getWidth());
                m_Locator.m_View.syncToLocator();
                el.applyStyles("width:0px;display:none;");
                this.pingClick('1-0-1');
            }
        }
        e.stopPropagation();
    },

    onClick: function(e) {
        var el =  Ext.get(this.m_selectionDivId);
        el.setWidth(0);
        el.hide();
        this.m_ResizeAction = false;
        
    },

    onContextMenu: function(e) {
        var menu = new Ext.menu.Menu();
        if (this.m_deferredContext) {
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
            this.m_ContextMenu = menu;
        } else {
            this.m_XY = e.getXY();
            e.preventDefault();  // this prevents the default contextmenu to open in Firefox (linux)
            e.stopPropagation();
        }
        
        this.pingClick('1-1');      
        var x_pos = this.m_XY[0] - Ext.get(this.m_DivId).getX()
        menu.add({text:'Set New Marker At Position', iconCls:'xsv-markers', scope:this,  
                handler:function() { this.m_App.newMarkerDlg(this, x_pos); this.pingClick('1-1-1');} 
        });
        menu.add('-');
        if (this.m_App.m_Origin) {
            menu.add({text:'Reset Sequence Origin', iconCls:'xsv-origin', scope:this, 
                handler:function() { this.m_App.clearOrigin(); this.pingClick('1-1-2'); }
            });
        } else {
            menu.add({text:'Set Sequence Origin At Position', iconCls:'xsv-origin', scope:this, 
                handler:function() { this.m_App.setOrigin(this, x_pos); this.pingClick('1-1-2'); }
            });
        }
        var nvMarker = this.m_App.getMarkersInfo().findMarkerByName(SeqView.MarkerNav);
        if (nvMarker) {
            menu.add({text:'Remove Navigation Marker', iconCls: 'xsv-marker-remove', scope:this,  
                handler:function() { nvMarker.deleteMarker(); this.pingClick('1-1-3'); }
            });
        }
        if (menu.items.length>0) { menu.showAt(this.m_XY); }
    },
    
    onResize: function(e) {},
    
    createMarkerElem: function(marker) {
        var elem = Ext.get(this.m_DivId);
        var create_params = marker.getCreateParams(true,this.m_Idx);
        var marker_elem = create_params.template.append(elem, create_params.options, true);
        marker_elem.setTop(this.m_TopOffset + 3);
        return marker_elem;
    }
    

});

//////////////////////////////////////////////////////////////////////////
SeqView.PanoramaSelection = function(panoramaDivId, selectionDivId, range) {
    this.m_selectionDivId = selectionDivId;
    this.m_panoramaDivId = panoramaDivId;
    this.m_Resizing = false;
    this.range = range;
    this.element = Ext.get(this.m_selectionDivId);
};

SeqView.PanoramaSelection.prototype = {

    pageToViewX: function(x) {
        var div_xy = Ext.get(this.m_panoramaDivId).getXY();
        return x - div_xy[0];
    },

    resize: function(xSelection1, xSelection2) {
        var x1 = this.pageToViewX(xSelection1);
        var x2 = this.pageToViewX(xSelection2);
        var el = this.element;
        el.setLeft(x1 <= x2 ? x1 : x2);
        el.setWidth(Math.abs(x2 - x1));
    }
  
};
/*  $Id: textview.js 37053 2016-11-30 22:10:41Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko
 * File Description: Sequence View Panel
 */

 /********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.TextView
/********************************************************************/


SeqView.TextView = (function() {
    var sm_TextViewRows = 100; // Max of 100 rows in the sequence text panels
    var sm_ViewTmpl = new Ext.Template(
        '<table width="{width}" class="xsv-fasta-table sv-fasta-table" style="border-style:solid;border-width:{border}px;" >',
        '<tbody><tr bordercolor="#ffffff"><td id="left_sequence_{id}" width=50 valign="top" align="right"><br><NOBR><em>{left_nums}</em></NOBR></td>',
        '<td valign="top"><span id="top_sequence_{id}"><em>{top_nums}</em></span><br/><div id="text_sequence_{id}">{sequence}</div></td>',
        '</tr></tbody></table>'
    );

    function constructor(app) {
        var div_id = 'seq-dlg_' + app.m_Idx;
        if (Ext.get(div_id)) { return; }
        this.m_App = app;
        this.m_Idx = 'S' + this.m_App.m_Idx; // m_Idx is needed for Locator object;
        //this.m_View = view;
        this.m_DivId = div_id;
        this.m_Color = '000000';
        this.m_FromSeq = 0;
        this.m_LenSeq = 0;
        this.m_ShowTranslation = true;
        this.m_SequenceCols = Ext.isWindows ? 90 : 100; // 10 rows, 10 bases in each row. 100 bases per line (configurable)
        this.m_SeqStarts = [];
        this.m_Locator = null;
        this.m_PrevXY = [];
        this.m_Flip = app.m_Flip;
        this.m_TextViewSize = this.m_SequenceCols * sm_TextViewRows; // load sequence size (for text views)
        this.m_SelectedFeat = '';
        this.m_HideTrans = 1;

    };

    constructor.getViewTmpl = function() {
        return sm_ViewTmpl;
    };
    constructor.getTextViewRows = function() {
        return sm_TextViewRows;
    };

    return constructor;

})();

////////////////////////////////////////////////////////////////
/// Public methods

SeqView.TextView.prototype = {

    canFlip: function() { return this.m_App.m_ViewParams.acc_type == 'DNA'; },
    isLoading: function() { return this.m_Loading; },
    isPanorama: function() { return false; },
    getURL: function() { return "&seq=" + (this.m_FromSeq+1) + ':' + (this.m_FromSeq+this.m_LenSeq) },
    toSeq: function() { return [this.m_FromSeq, this.m_FromSeq + this.m_LenSeq - 1, this.m_LenSeq]; },


    getTitle: function() {
        var title = "Sequence View";
        if (this.canFlip())
            title += (this.m_Flip ? " (negative" : " (positive") + " strand)";
        return title;
    },

    openDlg: function(obj) {
        if (this.m_SeqDlg) {
            this.m_SeqDlg.expand();
            return;
        }
        var from_seq, len_seq;
        // check if pass a marker position
        if (typeof obj == 'number') {
            from_seq = Math.max(0, obj - 949);
        } else if (Ext.isArray(obj)) {
            // range is passed from the URL
            from_seq = obj[0];
            len_seq = obj[1];
            this.m_Flip = obj[2];
        } else {
            if (!obj.isPanorama) {
                this.m_Flip = obj.orientation == 'reverse';
                this.displaySequenceTextDlg(obj.from, obj.to - obj.from, obj);
                return;
            }
            if (obj.isPanorama()) { // no parent
                from_seq = Math.max(0, this.m_App.m_SeqLength / 2 - this.m_TextViewSize / 2);
            } else {
                this.m_Flip = obj.getFlip();
                var vis_range = obj.toSeq();
                from_seq = Math.max(0, vis_range[0] + vis_range[2] / 2 - this.m_TextViewSize / 2);
            }
            from_seq = Math.round(from_seq/10) * 10;
        }
        len_seq = Math.min(from_seq + this.m_TextViewSize, this.m_App.m_SeqLength) - from_seq;
        this.displaySequenceTextDlg(from_seq, len_seq);
    },


    displaySequenceTextDlg: function(from_seq, len_seq) {
        this.m_App.resizeIFrame(500);
        if (typeof from_seq == 'number') {
            this.startImageLoading(from_seq, len_seq);
            if (this.m_App.m_Views[0].isPanorama()) {
                this.m_Locator = new SeqView.Locator(this, this.m_Color, false);
            }
        }
        this.m_SeqDlg = new Ext.Window({
            layout:'fit', title: this.getTitle(),
            app: this.m_App,
            width: 750, height: 450,  minWidth: 480,
            minHeight: 200,
            constrain:false,
            collapsible:true,
            maximizable:true,
            hideAction:'close',
            monitorResize: true,
            cls: 'SeqViewerApp',
            onEsc: function() { this.close(); },
            listeners:{
                'resize': function(window, width, height) {
                    var top_id = 'top_sequence_' + this.m_App.m_Idx;
                    var left_id = 'left_sequence_' + this.m_App.m_Idx;
                    var cur_top = document.getElementById(top_id);
                    if (!cur_top) { return; }
                    var starts_width = Ext.get(left_id).getWidth() + 40;
                    var line_width = cur_top.offsetWidth;
                    var one_width = line_width / this.m_SequenceCols;
                    this.m_SequenceCols = Math.floor((width - starts_width) / one_width / 10) * 10;
                    this.m_TextViewSize = this.m_SequenceCols * SeqView.TextView.getTextViewRows();
                    this.startImageLoading(this.m_FromSeq, Math.min(this.m_TextViewSize, this.m_App.m_SeqLength));
                },
                'maximize': function( window ){
                    window.setWidth(
                        Ext.isIE ? document.documentElement.clientWidth
                        : !Ext.isOpera ? document.body.clientWidth : self.innerWidth
                    );
                },
                'move': function(win, x, y) {
                    var scrollY = Ext.getDoc().getScroll().top;
                    var p = this.m_SeqDlg.el.translatePoints(x, scrollY);
                    if (y < p.top) {
                        this.m_SeqDlg.setPosition(x, p.top);
                    }
                },
                'beforeclose': function() { if (this.m_SeqDlg.maximized) this.m_SeqDlg.restore(); },
                'close': function() {
                    if (this.m_Locator) { this.m_Locator.remove(); }// remove locator form panoram
                    this.m_Locator = null;
                    this.m_App.un('origin_changed', this.onOriginChanged, this );
                    this.m_App.un('strand_changed', this.flipSeqPanelStrand, this );
                    this.m_App.m_TextView = null;
                    this.m_App.reCreateReflections();
                    delete this.m_SeqDlg;
                },
                scope: this
            },
            id: this.m_DivId,
            cls: 'SeqViewerApp',
            tbar:[
                {text:'Prev Page', scope: this, iconCls:'xsv-prev', handler:this.gotoSequencePrev},
                '-',
                {text:'Next Page', scope: this, iconCls:'xsv-next', handler:this.gotoSequenceNext},
                '->',
                //{text:'Show Gap', iconCls:'xsv-show_gap', handler:SeqApp.sequenceShowGap},'-',
                {xtype: 'combo',
                  itemId: 'tripletCombo',
                  queryMode: 'local',
                  forceSelection: 'true',
                  hidden: true,
                  triggerAction: 'all',
                  store: new Ext.data.ArrayStore({ fields: ['object_id', 'name'] }),
                  displayField: 'name',
                  valueField: 'object_id',
                  width: 250,
                  listeners: { select: function(combo) {
                      this.m_SelectedFeat = combo.getValue();
                      this.startImageLoading(this.m_FromSeq, this.m_LenSeq);
                  }, scope:this }
                },
                ' ',
                {tooltip:'Printer-Friendly Page', iconCls:'xsv-printer', scope:this, handler:function() { SeqView.PrinterFriendlyPage(this.m_App);/*SeqView.App.showPrintPageDlg(this.m_App.m_Idx);*/ }},
                '-',
                {text:'Flip Strands', id:'seq-flip-button', tooltip:'Flip Sequence Strands', iconCls:'xsv-flip-strands',
                    pressed:this.m_Flip, enableToggle:true, hidden: !this.canFlip(),
                    scope: this, handler:function() { this.m_App.setFlip(!this.m_Flip);/*this.flipSeqPanelStrand();*/} },
                this.m_App.m_ViewParams.acc_type != 'DNA' ? {hidden:true} : '-',
                {text: this.translationName(),
                   id: 'seq_trans_btn', //xtype: 'tbsplit',
                   tooltip:'Flip Translation', iconCls:'xsv-show_translation',
                   hidden: this.m_App.m_ViewParams.acc_type != 'DNA',
                   menuAlign:'tr-br?',
                   menu:new Ext.menu.Menu({items:[
                     {text:'Annotated', checked: this.m_HideTrans === 1, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 1; this.flipSeqTranslation(); } },
                     {text:'Conceptual', checked: this.m_HideTrans === 2, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 2; this.flipSeqTranslation(); } },
                     {text:'Both', checked: this.m_HideTrans === 0, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 0; this.flipSeqTranslation(); } },
                     {text:'None', checked: this.m_HideTrans === 3, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 3; this.flipSeqTranslation(); } }
                    ]})
                  },

                  this.m_App.m_ViewParams.acc_type != 'DNA'? {hidden:true} : '-',
                  {text:'Go To Position', iconCls:'xsv-goto_position', scope: this, handler:this.gotoSeqPositionDlg}

            ],
            modal: false,
            hidden: true,
            items:[{ xtype:'panel', autoScroll: true, id:'seq-panel_'+this.m_App.m_Idx, layout:'fit' }],
            buttons: [{text: 'Close', scope:this, handler: function() { this.m_SeqDlg.close(); }}]
        });
        var app = this.m_App;
        app.on('origin_changed', this.onOriginChanged, this );
        app.on('strand_changed', this.flipSeqPanelStrand, this );
        
        this.m_SeqDlg.on('close', function(){ app.resizeIFrame(); app.m_DialogShown = false; });

        app.m_DialogShown = true;
// This delay is not needed and in fact is buggy - when the callback for startImageLoading is called with the data
// the dialog can be not created yet.
//        if (app.m_iFrame) this.m_SeqDlg.show.defer(500, this.m_SeqDlg); else this.m_SeqDlg.show();
        this.m_SeqDlg.show();
    },

    onOriginChanged: function(app) {
        this.startImageLoading(this.m_FromSeq, this.m_LenSeq);
    },

    startImageLoading: function(from, len) {
        this.m_FromSeq = from;
        this.m_LenSeq = len;
        SeqView.App.simpleAjaxRequest({url: this.getSeqTextURL(), context:this, success: this.seqLoadCallback/*, error: this.seqLoadCallback*/});
    },

    getSeqTextURL: function() {
        var to_seq = this.m_FromSeq + this.m_LenSeq;
        var user_data_key = "";
        if (this.m_App.m_Key) {
            user_data_key = "&key=" + this.m_App.m_Key;
        }
        return this.m_App.m_CGIs.Sequence + '?opts=seq&id=' + this.m_App.GI + '&from=' + this.m_FromSeq + '&to=' + to_seq + '&col=' +
            this.m_SequenceCols + '&translation=' + (this.m_ShowTranslation ? 'true' : 'false') +
            '&reverse=' + (this.m_Flip ? 'true' : 'false') + '&selected=' + this.m_SelectedFeat + user_data_key;
    },

    genTopNumbers: function() {
        var top_nums = "";
        for (var i = 0; i != this.m_SequenceCols / 10; i++) {
            if (this.m_Flip) { top_nums += '<span>0</span>987654321'; }
            else { top_nums += '123456789<span>0</span>'; }
        }
        return top_nums;
    },

    seqLoadCallback: function(data, text) {
        var combo = this.m_SeqDlg.down('#tripletCombo');
        if (!combo) return;
        var from_cgi = data;

        var starts = from_cgi.starts
        var fromindex = 0;
        this.m_SeqStarts = [];
        while (true) { // build vertical positions array
            var newindex = starts.indexOf("<", fromindex);
            if (newindex == -1) { break; }
            var a_char = starts.substr(newindex + 1, 1);
            if (a_char=="d") {
                var tr_type = starts.substr(newindex+16,13);
                //! HACK
                if( this.m_SeqStarts.length > 0 ){
                    if (tr_type == 'seqtrans_prot') {
                        this.m_SeqStarts.push(-2);
                    } else {
                        this.m_SeqStarts.push(-1);
                    }
                }
                fromindex = newindex + 6;
            } else if (a_char=="b") {
                var num = parseInt( starts.substr(fromindex, newindex - fromindex) );
                this.m_SeqStarts.push(num);
                fromindex = newindex + 4;
            } else if (a_char=="/") {
                fromindex = newindex + 6;
            } else {
                break;
            }
        }

        if (this.m_App.getOrigin() != 0  ||  this.m_App.getFlip()) {
            var new_starts = '';
            re = /([^0-9]*)([0-9]+)(<br>)/ig;
            var arr;
            while ((arr = re.exec(starts)) != null) {
                new_starts += (arr[1] + this.m_App.posToLocal(arr[2]) + arr[3]);
            }

            starts = new_starts;
        }
        var top_nums = this.genTopNumbers();

        var html = SeqView.TextView.getViewTmpl().apply({border:0, width:'95%', left_nums:starts, top_nums:top_nums, id: this.m_App.m_Idx, sequence:from_cgi.sequence});
        Ext.getCmp('seq-panel_'+this.m_App.m_Idx).update(html);
        // Add tooltips
        var idx = 1;
        var ttl_search = from_cgi.app_id? 'ttl-'+from_cgi.app_id+'-' : 'ttl-id';
        var app = this.m_App;
        while (true) {
            var tip_elem = Ext.get(ttl_search + idx);
            if (!tip_elem) { break; }
            var idx_array = tip_elem.dom.getAttribute('rel').split(',');
            var feat_idx = idx_array[idx_array.length-1];
            var sig = "";
            for (j=0; j!=from_cgi.features.length; j++) {
                var feat = from_cgi.features[j];
                if (feat.id == feat_idx) { sig = feat.object_id; break; }
            }
            if (sig.length > 0) {
                var qtip = new Ext.ToolTip({
                    target: tip_elem,
                    dismissDelay: 5000,
                    cfg: { url: app.m_CGIs.ObjInfo+'?objinfo=1&signatures='+sig+'&id='+app.GI,
                        success: function(data) { this.setHtml(data[0].text); }
                    }
                });
                qtip.cfg.context = qtip;
                qtip.on('render', function(){ SeqView.App.simpleAjaxRequest(this.cfg); }, qtip);
            }
            idx++;
        } // while

        this.m_App.reCreateReflections(); // show reflections

        // add markers
        this.m_App.forEachMarker( function(m) { this.AddOrUpdateMarker(m); }, this);

        this.flipSeqTranslation();

        if (from_cgi.triplets) {
            combo.getStore().loadData(from_cgi.triplets.triplets);
            this.m_SelectedFeat = from_cgi.triplets.triplets[from_cgi.triplets.selected || 0].object_id;
            combo.setValue(this.m_SelectedFeat);
            combo.show();
            combo.collapse();
        } else {
            combo.hide();
            this.m_SelectedFeat = '';
        }

        Ext.get(this.m_DivId).on({
            'mousedown':    this.onMouseDown,
            'mouseup':      this.onMouseUp,
            'mousemove':    this.onMouseMove,
            'contextmenu':  this.onContextMenu,
            scope: this
        });

        this.m_App.updateLocator(this);
    },


    AddOrUpdateMarker: function(the_marker) {//num, color, seq_pos, name) {
        var seq_elem = Ext.get('text_sequence_' + this.m_App.m_Idx);
        if (seq_elem) {
            var char_height = this.getLineHeight();
            var char_width = this.getCharWidth();
            if (the_marker.span) {
                var marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+the_marker.marker_num+'_b');
                if (!marker) {
                    var tmpl = new Ext.Template('<div orig_id="marker_0_{num}" id="seqtext-marker_{aidx}_{num}_b" data-qtip="start of {name}" class="xsv-text_marker" style="border-color:{color};border-{side}: none; height:{height}; width:{width};"><div data-qtip="Locked" id="seqtext-marker_{aidx}_{num}_lock_b" style="display:{lock};" class="marker_locked_ov"></div> </div>');
                    marker = tmpl.append(seq_elem,
                                    {num:the_marker.marker_num, aidx:this.m_App.m_Idx,
                                     color:the_marker.color, name:the_marker.marker_name,
                                     lock:the_marker.lock?'inline':'none', side:this.m_App.getFlip()?'left':'right'}, true);
                }
                this.setTextMarker(marker, the_marker.seq_pos);
                marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+the_marker.marker_num+'_e');
                if (!marker) {
                    var tmpl = new Ext.Template('<div orig_id="marker_0_{num}" id="seqtext-marker_{aidx}_{num}_e" data-qtip="end of {name}" class="xsv-text_marker" style="border-color:{color};border-{side}: none;height:{height}; width:{width};"><div data-qtip="Locked" id="seqtext-marker_{aidx}_{num}_lock_e" style="display:{lock};" class="marker_locked_ov"></div> </div>');
                    marker = tmpl.append(seq_elem,
                                    {num:the_marker.marker_num, aidx:this.m_App.m_Idx,
                                     color:the_marker.color, name:the_marker.marker_name,
                                     lock:the_marker.lock?'inline':'none', side:this.m_App.getFlip()?'right':'left'}, true);
                }
                this.setTextMarker(marker, the_marker.seq_pos + the_marker.span);
            } else {
                var marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+the_marker.marker_num);
                if (!marker) {
                    var tmpl = new Ext.Template('<div orig_id="marker_0_{num}" id="seqtext-marker_{aidx}_{num}" data-qtip="{name}" class="xsv-text_marker" style="border-color:{color};height:{height}; width:{width};"><div data-qtip="Locked" id="seqtext-marker_{aidx}_{num}_lock" style="display:{lock};" class="marker_locked_ov"></div> </div>');
                    marker = tmpl.append(seq_elem,
                                    {num:the_marker.marker_num, aidx:this.m_App.m_Idx,
                                     color:the_marker.color, name:the_marker.marker_name,
                                     lock:the_marker.lock?'inline':'none'}, true);
                }
                this.setTextMarker(marker, the_marker.seq_pos);
            }
        }
    },


    UpdateLockMarker: function(marker) {
        if (marker.span) {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_lock_b');
            if (seq_marker) {
                seq_marker.setStyle('display', marker.lock ? 'inline' : 'none');
            }
            seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_lock_e');
            if (seq_marker) {
                seq_marker.setStyle('display', marker.lock ? 'inline' : 'none');
            }
        } else {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_lock');
            if (seq_marker) {
                seq_marker.setStyle('display', marker.lock ? 'inline' : 'none');
            }
        }
    },


    RemoveMarker: function(marker) {
        if (marker.span) {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
            if (seq_marker)
                seq_marker.remove();
            seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
            if (seq_marker)
                seq_marker.remove();
        } else {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num);
            if (seq_marker)
                seq_marker.remove();
        }
    },


    SetSeqPos: function(marker) {
        if (marker.span) {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num);
            if (seq_marker)
                seq_marker.remove();
        } else {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
            if (seq_marker)
                seq_marker.remove();
            seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
            if (seq_marker)
                seq_marker.remove();
        }
        this.AddOrUpdateMarker(marker);
    },


    getLineHeight: function() {
        this.txtMetrix = this.txtMetrix || new Ext.util.TextMetrics('top_sequence_' + this.m_App.m_Idx);
        return this.txtMetrix.getHeight("A");
    },
    getCharWidth: function() {
        this.txtMetrix = this.txtMetrix || new Ext.util.TextMetrics('top_sequence_' + this.m_App.m_Idx);
        return this.txtMetrix.getWidth("A");
    },


    setTextMarker: function(marker_elem, seq_pos) {
    
        if (seq_pos < this.m_FromSeq  ||  seq_pos >= this.m_FromSeq + this.m_LenSeq - 1) {
            marker_elem.setVisible(false);
            return;
        } else marker_elem.setVisible(true);

        var starts_width = Ext.get('left_sequence_' + this.m_App.m_Idx).getWidth() + 2;
        var line_width = document.getElementById('top_sequence_' + this.m_App.m_Idx).offsetWidth;
        var one_width = line_width / this.m_SequenceCols;

        var one_height_seq = this.getLineHeight();
        var one_height_trans = one_height_seq + 2;

        // 0-based to 1-based position shift
        seq_pos += 1;
        
        var marker_y = 2;
        var row_start = 0;
        var i;
        if (this.m_Flip) {
            for (i = 0; i != this.m_SeqStarts.length; i++) {
                row_start = this.m_SeqStarts[i];
                if( row_start < 0 ){
                    if( (row_start != -this.m_HideTrans || this.m_HideTrans == 0) && this.m_HideTrans != 3 ){
                        marker_y += one_height_trans;
                    }
                } else {
                    marker_y += one_height_seq;
                }

                if( row_start < 0 ) continue;
                if( seq_pos >= row_start - this.m_SequenceCols  &&  seq_pos < row_start ) break;
            }
            pos_h = row_start - seq_pos;
            var left_pos = starts_width + one_width * pos_h - 1;
            if( Ext.isIE7 ) left_pos += 2;
            marker_elem.setLeft( left_pos );
        } else {
            for (i = 0; i != this.m_SeqStarts.length; i++) {
                row_start = this.m_SeqStarts[i];
                if (row_start < 0) {
                    if ((row_start != -this.m_HideTrans || this.m_HideTrans === 0) && this.m_HideTrans !== 3) {
                        marker_y += one_height_trans;
                    }
                } else {
                    marker_y += one_height_seq;
                }

                if (row_start < 0) { continue; }
                if (seq_pos >= row_start  &&  seq_pos < row_start + this.m_SequenceCols) { break; }
            }
            var left_pos = starts_width + one_width * (seq_pos - row_start)  - 1;
            if (Ext.isIE) left_pos += 2;
            marker_elem.setLeft(left_pos);
        }
        marker_elem.setTop(marker_y - 2);
    },


    SeqText2SeqPos: function(new_top, new_left) { // sequence 2 screen (text view)
        var starts_width = Ext.get('left_sequence_' + this.m_App.m_Idx).getWidth() + 2;
        var line_width = document.getElementById('top_sequence_' + this.m_App.m_Idx).offsetWidth;
        var one_width = line_width / this.m_SequenceCols;

        var one_height_seq = this.getLineHeight();
        var one_height_trans = one_height_seq + 1;

        var seq_pos = 0;
        var x = 0;
        var y = 0;
        var row_start = 0;
        new_top = new_top - (one_height_seq + 2);
        var prev_row_start = -1;

        var half = one_height_seq / 2;

        for (var i = 0; i != this.m_SeqStarts.length; i++) {
            row_start = this.m_SeqStarts[i];
            if (row_start < 0) {
                if ((row_start != -this.m_HideTrans || this.m_HideTrans == 0) && this.m_HideTrans != 3) {
                    y += one_height_trans;
                }
                if (new_top <= y-half) { break; }
                else { continue; }
            }
            prev_row_start = row_start;
            if (new_top >= y-half  &&  new_top < y + one_height_seq-half) { break; }
            y += one_height_seq;
        }
        x = Math.ceil((new_left - starts_width) / one_width);

        if (this.m_Flip) { return prev_row_start - x - 1; }
        else { return prev_row_start + x; }
    },


    onMouseDown: function(e) {
        if (e.button === 0) {
            var elem = e.getTarget().id;
            if (elem.indexOf('seqtext-marker_') != 0) { return; }

            // check for lock
            var idx = elem.split('_')[2];
            var the_marker = this.m_App.findMarker('marker_0_' + idx);
            if (!the_marker || the_marker.lock) { return; }

            // ok, ready to drag
            this.m_PrevXY = e.getXY();
            this.m_CurMarker = Ext.get(elem);
            this.m_CurMarker.setStyle('cursor', 'move');
        }
    },


    onMouseUp: function(e) {
        if (!this.m_PrevXY || !this.m_CurMarker) { return; }

        var orig = this.m_CurMarker.dom.getAttribute('orig_id');
        var marker = this.m_App.findMarker(orig);

        var seq_pos = this.SeqText2SeqPos(this.m_CurMarker.getTop(true), this.m_CurMarker.getLeft(true));

        if (marker) {
            if (marker.span) {
                var eid = this.m_CurMarker.id;
                if (eid.indexOf('_b') != -1) {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
                    if (oelem) {
                        this.setTextMarker(oelem, seq_pos + marker.span -1);
                    }
                    marker.setSeqPos( seq_pos/*-1*/, false );
                } else {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
                    var pos = seq_pos - marker.span +1;
                    if (oelem) {
                        this.setTextMarker(oelem, pos);
                    }
                    marker.setSeqPos( pos, false );
                }
            } else {
                marker.setSeqPos( seq_pos/*-1*/, false );
            }
            this.setTextMarker(this.m_CurMarker, seq_pos);
        }

        this.m_PrevXY = null;
        this.m_CurMarker = null;
    },


    onContextMenu: function(e) {
        e.stopEvent();

        var elem = e.getTarget().id;
        //console.log(elem);
        if (elem.indexOf('seqtext-marker') == 0) {
            var menu = new Ext.menu.Menu();
            var idx = elem.split('_')[2];
            var marker = this.m_App.findMarker('marker_0_' + idx);
            if (marker) {
                menu.add(marker.marker_name);
                menu.add('-');
                menu.add({text:'Set To Position...', scope:this, handler:function() {
                    var cur_pos = marker.seq_pos-this.m_App.getOrigin();
                    if (cur_pos >= 0) { cur_pos++; }
                    Ext.MessageBox.prompt('Marker', 'Please enter new sequence position:', function(btn, text) {
                        if (btn!='ok' || text.length == 0) { return; }
                        var seq_pos = this.m_App.convertRelativePosition(text);
                        if (isNaN(seq_pos) || seq_pos < 0 || seq_pos >= this.m_App.m_SeqLength) {
                            Ext.MessageBox.show({title:m.marker_name, msg:'Invalid sequence position.', buttons:Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                            return;
                        }
                        marker.setSeqPos(seq_pos,true);
                    }, this, false, cur_pos ); }
                });
                menu.add('-');
                menu.add({iconCls:'xsv-marker-remove', scope:this, text:'Remove Marker', handler:function() { marker.deleteMarker(); } });
                menu.add({iconCls:'xsv-markers', text:'Marker Details', scope:this, handler:function() { this.m_App.showMarkersDlg( this ); } });

                menu.showAt(e.getXY());
            }
        }
    },


    onMouseMove: function(e) {
        if (!this.m_PrevXY || !this.m_CurMarker) {
            return;
        }

        var delta_x = this.m_PrevXY[0] - e.getXY()[0];
        var delta_y = this.m_PrevXY[1] - e.getXY()[1];
        this.m_PrevXY = e.getXY(); // save new values

        NCBIGBUtils.ClearBrowserSelection();

        var starts_width = Ext.get('left_sequence_' + this.m_App.m_Idx).getWidth() + 2;
        var line_width = document.getElementById('top_sequence_' + this.m_App.m_Idx).offsetWidth;

        var new_left = this.m_CurMarker.getLeft(true) - delta_x;
        var new_top = this.m_CurMarker.getTop(true) - delta_y;

        new_left = Math.min(new_left, starts_width + line_width-7);
        new_left = Math.max(new_left, starts_width);

        new_top = Math.min(new_top, Ext.get('text_sequence_' + this.m_App.m_Idx).getHeight());
        new_top = Math.max(new_top, Ext.get('top_sequence_' + this.m_App.m_Idx).getHeight()+1);

        var seq_pos = this.SeqText2SeqPos(new_top, new_left);

        this.m_CurMarker.setLeft(new_left);
        this.m_CurMarker.setTop(new_top);

        var orig = this.m_CurMarker.dom.getAttribute('orig_id');
        var marker = this.m_App.findMarker(orig);

        if (marker) {
            if (marker.span) {
                var eid = this.m_CurMarker.id;
                if (eid.indexOf('_b') != -1 ) {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
                    if (oelem) {
                        this.setTextMarker(oelem, seq_pos + marker.span -1);
                    }
                    marker.setSeqPos( seq_pos/*-1*/, false );
                } else {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
                    var pos = seq_pos - marker.span +1;
                    if (oelem) {
                        this.setTextMarker(oelem, pos);
                    }
                    marker.setSeqPos( pos, false );
                }
            } else {
                marker.setSeqPos( seq_pos/*-1*/, false );
            }
        }
    },


    gotoSequenceNext: function() {
        this.stepSequence(!this.m_Flip);
    },


    gotoSequencePrev: function() {
        this.stepSequence(this.m_Flip);
    },


    stepSequence: function(dir_fwd) {
        var step = dir_fwd ? this.m_TextViewSize : -this.m_TextViewSize;
        var from_seq = Math.min(this.m_FromSeq + step, this.m_App.m_SeqLength - this.m_TextViewSize);
        from_seq = Math.max(0, from_seq);
        this.startImageLoading(from_seq, this.m_LenSeq);
    },


    translationName: function() {
        if (this.m_HideTrans == 1)
            return "Annotated";
        else if (this.m_HideTrans == 2)
            return "Conceptual";
        else if (this.m_HideTrans == 3)
            return "None";
        else
            return "Both";
    },


    flipSeqPanelStrand: function() {
        this.m_Flip = !this.m_Flip;
        this.startImageLoading(this.m_FromSeq, this.m_LenSeq);
        Ext.getCmp('seq-flip-button').toggle(this.m_Flip);
        if (this.m_SeqDlg)
            this.m_SeqDlg.setTitle(this.getTitle());
    },


    flipSeqTranslation: function() {
        var items, i, el;
        if (this.m_HideTrans == 0 || this.m_HideTrans == 3) {
            items = Ext.query('div[class^=xsv-seqtrans_]');
            for(i = 0; i < items.length; i++) {
                el = Ext.get(items[i]);
                el.setStyle('display',(this.m_HideTrans == 0 ? 'block' : 'none'));
                if (this.m_HideTrans == 0) {
                    Ext.util.CSS.updateRule('.xsv-ttu3','border-top-width', '0px');
                    Ext.util.CSS.updateRule('.xsv-ttu7', 'border-top-width', '0px');
                }
            }
        } else {
            Ext.util.CSS.updateRule('.xsv-ttu3', 'border-top-width', '2px');
            Ext.util.CSS.updateRule('.xsv-ttu7', 'border-top-width', '2px');
            items = Ext.query('div[class^=xsv-seqtrans_prot]');
            for(i = 0; i < items.length; i++) {
                el = Ext.get(items[i]);
                el.setStyle('display',(this.m_HideTrans == 1 ? 'block' : 'none'));
            }
            items = Ext.query('div[class^=xsv-seqtrans_trans]');
            for(i = 0; i < items.length; i++) {
                el = Ext.get(items[i]);
                el.setStyle('display',(this.m_HideTrans == 2 ? 'block' : 'none'));
            }
        }
        this.m_App.forEachMarker(function(m) {
            if (m.span) {
                var text_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+m.marker_num+'_b');
                this.setTextMarker(text_marker, m.seq_pos);
                text_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+m.marker_num+'_e');
                this.setTextMarker(text_marker, m.seq_pos + m.span - 1);
            } else {
                var text_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+m.marker_num);
                this.setTextMarker(text_marker, m.seq_pos);
            }
        },this);
        el = Ext.getCmp('seq_trans_btn');
        el.setText(this.translationName());
    },


    gotoSeqPositionDlg: function() {
        Ext.MessageBox.prompt('Go to position',
        'Please enter sequence position<br />(possible&nbsp;value&nbsp;formats&nbsp;are&nbsp;10k,&nbsp;-20,&nbsp;1m):',
            function(btn, text) {
            if (btn!='ok'  || text.length === 0) { return; }
            var position = 0;
            var bad_pos = !NCBIGBUtils.isNumeric(text);
            if (!bad_pos) {
                position = NCBIGBUtils.stringToNum(text);
                position = this.m_App.posToGlobal(position)
                if (position >= this.m_App.m_SeqLength) {
                    bad_pos = true;
                }
            }
            if (bad_pos) {
                Ext.MessageBox.alert('Go to position', 'Invalid position.');
            } else {
                var from_seq = Math.min(position, this.m_App.m_SeqLength - this.m_TextViewSize);
                from_seq = Math.max(0, from_seq);
                this.startImageLoading(from_seq, this.m_LenSeq);
            }
        }, this, false);
    },


    syncToLocator: function() {
        if (!this.m_Locator) {
            return;
        }
        var from_seq = this.m_App.m_Panorama.toSeq( this.m_Locator.getLeft(true) );
        this.startImageLoading(from_seq, this.m_LenSeq);
    }
};

/*  $Id: alignview.js 37945 2017-03-07 20:15:44Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko
 * File Description:
 */
 
 
/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.AlignSelection 
/********************************************************************/

SeqView.AlignSelection = function(view, area) {
   this.area = area;
   this.m_View = view; // view index
   //this.signature = area['signature'];
   this.parent_elem = Ext.get(this.m_View.m_DivId);
   
   var tpl = new Ext.Template('<div id="align_selection_id_{idx}" class="over_selection"/>');
   
   this.element = tpl.append(this.parent_elem, {idx:this.m_View.m_Idx}, true);
   this.element.setLeft(area.x); // the output from cgi seems a bit off. That's why this offsets. Or is it browser-specific? :)
   this.element.setTop(area.y); 
   this.element.setHeight(area.h);
   this.element.setWidth(area.w);
   
   this.element.on({
        //'dblclick' : SeqApp.onGViewDblClick,
        'click' :   this.m_View.onClick,
        //'mousedown' : SeqApp.onGViewMouseDown,
        //'mouseup' : SeqApp.onGViewMouseUp,
        'mousemove' : function(e) {e.stopPropagation();} ,
        //'contextmenu' : SeqApp.onGViewContextMenu
        scope: this.m_View
   });
   
    this.qtip = new Ext.ToolTip({
        target: this.element, 
        trackMouse:false, 
        autoWidth:true, 
        autoHide:true, 
        html:area.descr, 
        dismissDelay:5000,
        cls: 'SeqViewerApp'
    });
};


SeqView.AlignSelection.prototype = {
    movePix: function(delta) {
        new_left = this.element.getLeft(true) - delta;
        this.element.setLeft(new_left);
    },
    remove: function() { 
        if (this.qtip) {
            this.qtip.destroy();
        }
        this.element.remove(); 
    }
};


/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Alignment 
/********************************************************************/

Ext.define('SeqView.Alignment', { 
    extend: 'SeqView.View',
    m_PrevXY: null,
    m_FromSeq: 0,
    m_LenSeq:  0,
    m_ViewPageSize: 500,
    m_Expand: '',
    m_History: [],

    constructor: function(app) { this.callParent(['alignment', app]); },

    isAlignment: function() { return true; },

    createPanel: function() {
        this.m_DivId = 'alignment_id' + this.m_Idx;
        
        var tbar = [];
        
        if (this.m_App.m_Embedded === false) {
            this.m_Spacer = this.m_App.addView({
				border: false, 
				height: this.getSpacerHeight()
			}); // spacer
            tbar.push(this.createChooseColorBtn(), '-');
        }
        var range = this.toSeq();
        var r_range_from = this.m_App.posToLocal(range[0]);
        var r_range_to   = this.m_App.posToLocal(range[1]);
        if( this.m_App.m_Flip ){
            var t = r_range_from;
            r_range_from = r_range_to;
            r_range_to = t;
        }
        
        if( this.m_App.m_Toolbar["history"] == true ){
            tbar.push(
                { iconCls:'back', tooltip: 'Go back', itemId: 'histPrev'+this.m_Idx, scope:this, 
                    handler:function() { 
                        this.stepHistory(); 
                    }
                }
            );          
        }
        
        if( this.m_App.m_Toolbar["name"] == true ){
            var tiptitle = '<b>'
            tiptitle += this.m_App.m_ViewParams['id'];
            tiptitle += ': ' + this.m_App.m_Config.SeqInfo.title;
            tiptitle += '</b><br/><br/>';
            
            tiptitle += r_range_from.commify() + '&nbsp;-&nbsp;' + r_range_to.commify();
             
            tiptitle += '&nbsp;('+ (range[1] - range[0] + 1).commify() + '&nbsp;';
            tiptitle += (this.m_App.m_ViewParams['acc_type']=='protein' ? 'residues' : 'bases');
            tiptitle += '&nbsp;shown';
            
            if( this.m_App.m_ViewParams.acc_type == 'DNA' ){ 
                if( this.m_App.m_Flip ){
                    tiptitle += ',&nbsp;negative strand'; 
                } else { 
                    tiptitle += ',&nbsp;positive strand'; 
                }
            }
            tiptitle += ")";
            
            var title = '<b>' + r_range_from.commify() + ' - ' + r_range_to.commify() + ' (' + range[2].commify() + ' bases shown)</b>';
            
            
            tbar.push({ text: title, itemId: 'tbtitle'+this.m_Idx, tooltip: tiptitle, width: 200, handler: this.openFullView,  scope:this });
        }
        if( this.m_App.m_Toolbar["panning"] == true ){
                tbar.push(
                '-',
                    {iconCls:'pan-left', tooltip: 'Pan left', scope:this, handler:function() { this.gotoPrev(); }},
                    {iconCls:'pan-right', tooltip: 'Pan right', scope:this, handler:function() { this.gotoNext(); }},
                    '-'
                );
            }                       
                            
        if( this.m_App.m_Toolbar["zoom"] == true ){
            tbar.push(
                {text: '&nbsp;&nbsp;-&nbsp;&nbsp;', scale: 'small',tooltip:'Zoom Out', scope:this, handler:function() { this.zoomOut();}}
            );
        
            this.m_TbSlider = new Ext.Slider({
                name: 'zoom', width: 100,  
                minValue: 0, maxValue: 100, value: 0, increment: 5,
                topThumbZIndex: 100,
                tipText: function(thumb){ return String(thumb.value) + '%'; },
                listeners: { 
                    changecomplete: function(el,val) { 
                      var slider_range = 100;
                            var slider_val = 100 - val;
                                      
                            var vis_range = this.toSeq();
                            var zoom_point = vis_range[0] + vis_range[2] / 2;
                                      
                            var max_bases_to_show = this.m_App.m_SeqLength;
                            var rightPix = this.m_View.body.getWidth() - 2;
                            var min_bases_to_show = rightPix * SeqView.MinBpp;
                       
                            var ee = Math.log(max_bases_to_show / min_bases_to_show);
                             
                            var slider_val = Math.min(slider_range, Math.max(0, slider_val));
                             
                            var len = min_bases_to_show * Math.exp(ee * slider_val / slider_range);
                            var from = zoom_point - len / 2;

                            len = Math.min(this.m_App.m_SeqLength, Math.round(len));
                            from = Math.round( Math.max(0, from) );
                            if (from + len > this.m_App.m_SeqLength) {
                                var extra = from + len - this.m_App.m_SeqLength;
                                from -= extra;
                            }
                            this.pushToHistory(range[0], range[1]-range[0],"");
                            this.loadImage(from, len);

                    }.bind(this)
                }
            });
            tbar.push( this.m_TbSlider );   

            tbar.push(
                {text: '&nbsp;&nbsp;+&nbsp;&nbsp;', scale: 'small',tooltip:'Zoom In', scope:this, handler:function() { this.zoomIn();} },
                {iconCls:'xsv-zoom_seq', tooltip:'Zoom To Sequence', scope:this, handler:function() { this.zoomSeq();} }
            );
        }
        tbar.push('->');
        if( this.m_App.m_Toolbar["reload"] == true ){
            tbar.push('-',
                {iconCls:'x-tbar-loading', id:'seq-view-loading-'+this.m_Idx, tooltip:'Reload View', 
                    scope: this, handler:function(){this.refresh(); }, region: 'east'//, style: {"position": "relative", "right":"0"} 
                }
            );
        }
        
        var first_graphic = true;
        for( var i = 0; i < this.m_App.m_Views.length; i++ ){
            var view = this.m_App.m_Views[i];
            if( view && view != this && view.isGraphic() ){
                first_graphic = false;
                break;
            }
        }
        if( first_graphic ){
            if( this.m_App.m_Toolbar["help"] == true ){
            tbar.push(
                {
                    iconCls: 'xsv-question', 
                    tooltip:'Help', 
                    layout: {type:'vbox'},
                    scope: this, 
                    menu:new Ext.menu.Menu({
                        defaultOffsets: [-70,0],
                        scope: this,
                        items:[{
                            text: 'Help',
                            iconCls: 'xsv-question', 
                            scope:this,
                            handler: function() {SeqView.showHelpDlg(); this.pingClick('2-2-10-1');}
                        },{
                            text: 'Link to View', 
                            iconCls: 'xsv-link_to_page',
                            scope: this, 
                            handler: function() {this.m_App.showLinkURLDlg(); this.pingClick('2-2-10-3');}
                        },{
                            text: 'Feedback', 
                            iconCls: 'xsv-feedback',
                            scope: this, 
                            handler: function() {this.m_App.showFeedbackDlg(); this.pingClick('2-2-10-0');}
                        },{
                            text: 'About', 
                            scope: this, 
                            handler: function() {SeqView.showAboutMessage(this.m_Idx); this.pingClick('2-2-10-2');}
                        }]
                    })
                }
            );
            }
         } else {
            tbar.push(
                {iconCls: 'tb-close', tooltip:'Close', scope: this, handler:function() { this.remove(); }}
            );
        }
        
        var tools = [];
        if (this.m_App.m_Embedded === false) {
            tools.push(this.createCloseTBar());
        } else {
            tools.push({id:'help',qtip: 'Help', handler: function(){ SeqView.showHelpDlg(); }})
        }
        
        this.m_View = this.m_App.addView({
            //collapsible:true, title:'Loading...',
            collapsible:false, header: false,
            tools:tools, view:this, tbar:tbar, html:'<div class="alignment_div" id="'+this.m_DivId+'"/>'
        });
        
        
        Ext.get(this.m_DivId).on({
            'mousemove' : this.onMouseMove,
            'contextmenu':this.onContextMenu,
            scope: this
        });
        
        if (this.m_App.m_Embedded === false)
            this.m_Locator = new SeqView.Locator(this, this.m_Color, true);

        this.loadImage(this.m_FromSeq, this.m_LenSeq);
    },
    
    loadImage: function(from,len) {
        this.m_Loading = true; // start loading
        this.m_Width = this.getScreenWidth();
      
        var params = {id: this.m_App.GI, view:'aln', client:'assmviewer', width: this.m_Width, from: from, len:len}
        if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0)
            params.appname = this.m_App.m_AppName;
        var url = this.m_App.m_CGIs.Alignment;
        if (this.m_Expand && this.m_Expand.length>0) params.expand = this.m_Expand;
        this.addURLParams(params);

        // Add all 'data_key' keys from all 'alignment_track's to the 'key' parameter
        // to feed the data source
        var tracks = this.m_App.getTrackObjects();
        var key = params.key || "";
        for (var i = 0; i < tracks.length; ++i) {
            var track = tracks[i];
            if (track.key === 'alignment_track' && track.data_key) {
                if (key) key += ',';
                key += track.data_key;
            }
        }
        if (key) params.key = key;

      
        Ext.getCmp('seq-view-loading-'+this.m_Idx).disable();
        this.m_App.AjaxRequest({url: url, context: this, data: params,
            success:this.checkJobStatus, error:this.loadFailure});      
    }, 

    checkJobStatus: function(data, text, res) {
        if (!this.m_Loading) return;

        var from_cgi = SeqView.decode(data);
        if (from_cgi.job_status) {
            if (from_cgi.job_status == 'failed') {
                this.loadFailure(null, from_cgi.error_message);
            } else if(from_cgi.job_status == 'canceled') {
                this.loadFailure(null, 'Job canceled');
            } else {
                var url = this.m_App.m_CGIs.Alignment + '?job_key=' + from_cgi.job_id
                Ext.defer(this.m_App.AjaxRequest, 2000,this,[{url:url, context: this,
                        success: this.checkJobStatus, error: this.loadFailure}]);
            }
        } else {
            Ext.getCmp('seq-view-loading-'+this.m_Idx).enable();
            if (from_cgi.error) {
                this.loadFailure(null, from_cgi.error);
            } else if (from_cgi.success === false) {
                this.loadFailure(null, from_cgi.msg);
            } else {
                // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
                // This is a way to provide reliable URL resolution for embedding. SV-1760
                if (from_cgi.img_url && from_cgi.img_url.charAt(0) == '?') {
                    from_cgi.img_url = this.m_App.m_CGIs.NetCache + from_cgi.img_url;
                }
                this.m_FromCgi = from_cgi;
                this.m_FromSeq =      this.m_FromCgi.from; // got from cgi 0-based
                this.m_LenSeq =       this.m_FromCgi.len;
                this.m_AlgnLen =      this.m_FromCgi.seq_length;
                this.m_ViewPageSize = this.m_LenSeq;
                this.m_Width =        this.m_FromCgi.img_width;
                this.m_Height =       this.m_FromCgi.img_height;
                this.m_Loading =      false; // loaded
    
                var the_div = Ext.get(this.m_DivId);
                var img_el = the_div.first('img');
                if (img_el) {
                    var d = Ext.getDom(img_el);
                    d.src = from_cgi.img_url;
                } else {
                    var tpl = new Ext.Template('<img class="sv-drag sv-highlight sv-dblclick" src="{img_url}">'),
                    img_el = tpl.append(the_div, from_cgi,true);
                }
                the_div.setStyle('height', this.m_Height + 'px');
                this.m_View.updateLayout();
                if( this.m_TbSlider ){  
                    var slider_range = 100;

                    var range = this.toSeq();
                    this.pushToHistory(range[0], range[1]-range[0],"");
                    var vis_len = range[2];
                              
                    var max_bases_to_show = this.m_App.m_SeqLength;
                    var rightPix = this.m_View.body.getWidth() - 2;
                    var min_bases_to_show = rightPix * SeqView.MinBpp;
                    //var min_bases_to_show = this.m_App.getPanoramaWidth() / 10;
                              
                    var ee = Math.log(max_bases_to_show / min_bases_to_show);
                     
                    var slider_val = slider_range * Math.log(vis_len / min_bases_to_show) / ee;

                    this.m_TbSlider.setValue( 100 - slider_val );
                }
                this.m_App.notifyViewLoaded(this);
            }            
        }
    },
    loadFailure: function(data, res) {
        Ext.MessageBox.show({title:'Image loading error', msg:res, buttons:Ext.MessageBox.OK, icon:Ext.MessageBox.INFO});
        var the_div = Ext.get(this.m_DivId);
        the_div.setStyle('background-image', 'none');
        this.m_View.setTitle("Not available");
    },

    onMouseDown: function(e) {},
    onMouseUp: function(e) {},
    onMouseMove: function(e) {
        var area = this.hitTest(e.getXY());
        
        var elem = Ext.get( this.m_DivId );
        
        if (area && !this.m_Selection) {
            elem.setStyle('cursor', 'pointer');
            this.m_Selection = new SeqView.AlignSelection(this, area);
        } else if (this.m_Selection) {
            elem.setStyle('cursor', 'default');
            this.m_Selection.remove();
            this.m_Selection = null;
        }
    },

    onClick: function(e)
    {
        var area = this.hitTest(e.getXY());
        if (area) { // open alignment
            var cur_name = area.name;
            var cur_expand = this.m_Expand;
            var e_array = cur_expand.length > 0 ? cur_expand.split(',') : [];
            var expanded = false;

            for (var i = 0;  i != e_array.length;  i++) {
                if (e_array[i] == cur_name) {
                    expanded = true;
                    delete e_array[i];
                    break;
                }
            }
      
            if (!expanded) e_array.push(cur_name);
            this.m_Expand = e_array.join(',');      
            this.loadImage(this.m_FromSeq, this.m_LenSeq);
        }
    },

    onMouseOut: function(e) {
    },

    onContextMenu: function(e) {
        e.stopEvent();
        menu = new Ext.menu.Menu();
             
        menu.add({iconCls:'xsv-zoom_plus', text:'Zoom In', scope:this, handler:function() { this.zoomIn();} });
        menu.add({iconCls:'xsv-zoom_minus', text:'Zoom Out', scope:this, handler:function() { this.zoomOut();} });
        menu.add('-');
        menu.add({iconCls:'xsv-zoom_seq', text:'Zoom To Sequence', scope:this, handler:function() { this.zoomSeq();} });
        menu.add('-');
        menu.add({iconCls:'expand_all', text:'Expand All', scope:this, handler:function() { this.expandAll();} });
        menu.add({iconCls:'collapse_all', text:'Collapse All', scope:this, handler:function() { this.collapseAll();} });
        menu.add('-');
        menu.add({iconCls:'prev', text:'Prev Page', scope:this, handler:function() { this.gotoPrev();} });
        menu.add({iconCls:'next', text:'Next Page', scope:this, handler:function() { this.gotoNext();} });
        menu.showAt(e.getXY());
    },

    toSeq: function() { 
         return [this.m_FromSeq, this.m_FromSeq + this.m_LenSeq - 1, this.m_LenSeq];
    },
    
    updateTitle: function() {
        var range = this.toSeq();
        var r_range_from = range[0] + 1;
        var r_range_to   = range[1] + 1;

        var title = this.m_App.m_Embedded !== false ? "<a href='#' class='hidden_for_print' id='new_view_link_al"+this.m_App.m_Idx+"' style='float:right;padding-right:7px;'>Open Full View</a>" : 'Alignment View: ';
        title += r_range_from.commify() + ' - ' + r_range_to.commify() + ' (' + range[2].commify() + ' bases shown)';
        this.m_View.setTitle(title);

        var href = Ext.get('new_view_link_al'+this.m_App.m_Idx); 
        if (href) { 
            href.on({ 'click' : this.openFullView,  scope:this });
        }

        this.m_App.updateLocator(this);
    },

    hitTest: function(page_xy) {
        var elem_xy = Ext.get( this.m_DivId ).getXY();
        var xx = page_xy[0] - elem_xy[0];// - config['scroll_pix'];
        var yy = page_xy[1] - elem_xy[1];// - config['top_offset'];
        for (var a = 0; a < this.m_FromCgi.checkboxes.length; a++) {
            var area = this.m_FromCgi.checkboxes[a];
            var the_x = area.x;
            var the_y = area.y;
            var the_w = area.w;
            var the_h = area.h;
            if (xx >= the_x-1 && xx <= the_x+the_w && yy >= the_y-1 && yy <= the_y+the_h) return area; // add extra pixel
        }
        return null;
    },
    
    zoomIn: function() {
    	 var new_len  = this.m_LenSeq / 2;
        var new_from = this.m_FromSeq + new_len / 2;
        this.loadImage(new_from, new_len);
    },

    zoomSeq: function( center_seq_pos ){        
        
        if( !center_seq_pos ){
            center_seq_pos = Math.floor( this.m_FromSeq + this.m_LenSeq /2 );
        }

        var new_len  = Math.floor( this.getScreenWidth() * SeqView.MinBpp ); 
        var new_from = Math.floor( center_seq_pos - new_len /2 );
        if( new_from < 0 ){
            new_from = 0;
        }
        if( new_from + new_len > this.m_App.m_SeqLength ){
            new_len = this.m_App.m_SeqLength - new_from;
        }

        this.loadImage( new_from, new_len );
    },


    zoomOut: function() {
        var new_len  = this.m_LenSeq * 2;
        var new_from = this.m_FromSeq - this.m_LenSeq / 2;
   
        new_from = Math.max(0, new_from); // not to exteed the sequence range
        if (new_from + new_len > this.m_AlgnLen) new_len = this.m_AlgnLen - new_from;

        this.loadImage(new_from, new_len);
    },
    

    gotoPrev: function() {
        var new_from = Math.max(0, this.m_FromSeq - this.m_ViewPageSize);
        this.loadImage(new_from, this.m_LenSeq);
    },
    

    gotoNext: function() {
        var new_from = Math.min(this.m_FromSeq  + this.m_ViewPageSize - 1, this.m_AlgnLen - this.m_ViewPageSize);
        new_from = Math.max(0, new_from);
        this.loadImage(new_from, this.m_LenSeq);
    },  

    expandAll: function() {
        var e_array = [];
        for (var a = 0; a < this.m_FromCgi.checkboxes.length; a++) 
            e_array.push( this.m_FromCgi.checkboxes[a].name );
   
        this.m_Expand = e_array.join(',');   
        this.loadImage(this.m_FromSeq, this.m_LenSeq);
    },

    collapseAll: function() {
        this.m_Expand = '';
        this.loadImage(this.m_FromSeq, this.m_LenSeq);
    },
    
    refresh: function(options) {
        if (this.m_VisLenSeq == 0)
            return;
        this.loadImage(this.m_FromSeq,this.m_LenSeq);
    },


    syncToLocator: function() {
         if (!this.m_Locator) {
            return;
        }            
        var vis_from, vis_to;
        if (this.m_Locator.m_ResizeRight) { vis_from = this.m_FromSeq; }
        else { vis_from = this.m_App.m_Panorama.toSeq( this.m_Locator.getLeft(true) ); }

        if (this.m_Locator.m_Scroll) { vis_to = vis_from + this.m_LenSeq-1; }// keep length 
        else { vis_to = this.m_App.m_Panorama.toSeq( this.m_Locator.getLeft(true) + this.m_Locator.getWidth()+3 ); }
        this.loadImage(vis_from, vis_to-vis_from+1);
    },
    
    checkLocatorWidth: function(width) { return true; },
    
    openFullView: function() { this.m_App.getLinkToThisPageURL('portal'); },
    
    pushToHistory: function( from, length, title ){
        if (this.m_HistIx < 0) this.m_HistIx = 0;
        var current = { from: from, len: length, title: title };

        if (this.m_HistIx < this.m_History.length) {
            if (this.mf_HistoricalBack
                || (Math.abs(current.from - this.m_History[this.m_HistIx].from) <= 1
                    && current.len == this.m_History[this.m_HistIx].len)) {
                this.mf_HistoricalBack = false;
                return;
            }

            if( this.m_HistIx > 0 ){
                this.m_History = this.m_History.slice( this.m_HistIx );
                this.m_HistIx = 0;
            }
            
            this.m_History.unshift( current );
            
        } else {
            this.m_History.push( current );
            this.m_HistIx = this.m_History.length -1;
        }
        this.mf_HistoricalBack = false;
        
        var MAX_HISTORY_LENGTH = 10;
        var mhl = MAX_HISTORY_LENGTH;
        
        if( this.m_History.length > mhl ){
            if( this.m_HistIx < mhl ){
                this.m_History = this.m_History.slice( 0, mhl );
            } else {
                this.m_History = this.m_History.slice( this.m_HistIx -mhl +1, mhl );
                this.m_HistIx = this.m_History.length -1;
            }
        }
        this.updateHistoryButtons();
    },


    stepHistory: function( forward ){
        if( forward ){
            if( this.m_HistIx > 0 ){
                this.m_HistIx--;
                this.mf_HistoricalBack = true;
                this.loadImage( this.m_History[this.m_HistIx].from, this.m_History[this.m_HistIx].len ); 
            }
        } else {
            if( this.m_HistIx < this.m_History.length-1 ){
                this.m_HistIx++;
                this.mf_HistoricalBack = true;
                this.loadImage( this.m_History[this.m_HistIx].from, this.m_History[this.m_HistIx].len ); 
            }
        }
        this.updateHistoryButtons();
    },  


    updateHistoryButtons: function() {
        var tbar = this.m_View.getDockedComponent(0);
        if( tbar ){
            var btn_prev = tbar.getComponent( 'histPrev'+this.m_Idx );
            if( btn_prev ){
                var tooltip = '';
                if( this.m_HistIx >= this.m_History.length-1 ){
                    btn_prev.disable();
                    tooltip = 'Back';

                } else {
                    btn_prev.enable();
                    
                    if( this.m_HistIx+1 < this.m_History.length ){
                        tooltip += 'Back to ' + this.m_History[this.m_HistIx+1].title;
                    }
                }
                btn_prev.setTooltip( tooltip );
            }
            
            var btn_next = tbar.getComponent( 'histNext' );
            if( btn_next ){
                var tooltip = '';
                if( this.m_HistIx <= 0 ){
                    btn_next.disable();
                    tooltip = 'Forward';

                } else {
                    btn_next.enable();
                    
                    for( var j = this.m_HistIx-1; j >= 0; j-- ){
                        tooltip += this.m_History[j].title + '<br/>';
                    }
                }
                btn_next.setTooltip( tooltip );
            }
        }
    }

});

/*  $Id: graphview.js 38403 2017-05-04 19:38:17Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 * File Description:
 */


SeqView.ChunkWidth = 4000;
SeqView.MinBpp = 1/24;
SeqView.PreloadMargin = SeqView.ChunkWidth/10; // load next chunk when the current is with this.m_PreloadMargin to the edge

SeqView.PAN_RIGHT = 1;
SeqView.PAN_LEFT = -1;

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Graphic
/********************************************************************/

Ext.define('SeqView.Graphic', {
    extend: 'SeqView.View',
    m_bottomRulerHeight: 0,
    m_PrevCgi: null,
    m_NextCgi: null, // image caching
    m_FromSeq: 0,
    m_LenSeq:  0,
    m_VisFromSeq: 0,
    m_VisLenSeq: 0,
    m_ScrollPix: 0,
    m_Flip: false,
    m_Selection: null,
    m_RangeSelectionSet: [],
    m_Reflections: null,
    m_SelectedSig: null,
    m_SelectedRangeSet: [],
    m_TbSlider: null,
    m_DocMouseMove: null,
    m_DocMouseUp: null,
    m_DocTouchMove: null,
    m_DocTouchUp: null,

    constructor: function(app) {
        this.clear();
        SeqView.Graphic.superclass.constructor.apply(this, ['graphical',app]);
        this.m_LenSeq = this.m_App.m_SeqLength;
        this.m_Flip = app.getFlip();
    },

    isGraphic: function() { return true; },
    getSelectionTop: function() { return 0; }, // right from the top
    getSelectionHeight: function() { return this.m_Height; },// down to the bottom
    clearSelectedSig: function() { return this.setSelectedSig(""); },
    canFlip: function() { return this.m_App.m_ViewParams.acc_type == 'DNA'; },
    flipStrand: function() { this.setFlip(!this.getFlip()); },
    getFlip: function() { return this.m_Flip; },

    setFlip: function(flip) {
        if (flip != this.m_Flip) {
            this.setFlipNoReload(flip);
            this.refresh();
        }
    },

    setFlipLocal: function( flip ){
        if( flip != this.m_Flip ){
            this.m_Flip = flip;
            var button = Ext.getCmp( 'flip-button-' + this.m_Idx );
            if( button ) button.toggle( flip );
            this.refresh();
            this.pushToHistory();
        }
    },

    setFlipNoReload: function(flip,no_broadcast) {
        this.m_Flip = flip;
        var button = Ext.getCmp('flip-button-' + this.m_Idx);
        if (button) button.toggle(flip);
        if (typeof no_broadcast == "undefined" || !no_broadcast) {
            // Propagate to the peer views
            var this_view = this;
            this.m_App.setFlip(flip, this_view);
        }
    },

    getMostRightPix: function() { return this.getScreenWidth(); },

    seq2Pix: function(seq_pos) {
        var pos = 0;
        if (this.getFlip()) {
            pos = (this.m_LenSeq - (seq_pos - this.m_FromSeq) - 1 + 0.5) * this.m_bpPix;
        } else {
            pos = (seq_pos - this.m_FromSeq + 0.5) * this.m_bpPix;
        }
        return Math.round(pos);
    },

    seq2PixScrolled: function(seq_pos) { return this.seq2Pix(seq_pos) + this.m_ScrollPix; },

    pix2Seq: function(pix_pos, adj) { // screen 2 sequence (gview)
        if (adj === undefined) adj = 0;
        if (this.getFlip()) pix_pos = this.m_Width - pix_pos;   // FIXME
        return  this.m_FromSeq + Math.floor(pix_pos / this.m_bpPix + adj);

    },

    getScreenHeight: function() { return this.m_View.body.getHeight() - (this.m_App.m_Embedded ?  0 : 2); },
    getMinLength: function(){ return Math.floor(this.getScreenWidth() * SeqView.MinBpp); }, 

    toSeq: function() { // screen 2 sequence range (more precise calculation required!!!)
        if (!this.m_FromCgi)
            return [this.m_FromSeq, this.m_FromSeq + this.m_LenSeq - 1, this.m_LenSeq];

        // If we positioned exacly as requested and did not move,
        // use exact cached position
        if (this.m_ReqFrom !== undefined && this.m_ReqLen !== undefined)
            return [this.m_ReqFrom, this.m_ReqFrom + this.m_ReqLen - 1, this.m_ReqLen];

        var ppb = this.m_Width / this.m_LenSeq;
        var view_seq_width = Math.floor(this.getScreenWidth()/ppb);

        var from, to;
        if (this.getFlip()){
            var left = this.m_FromSeq + this.m_LenSeq - 1;
            to = left - Math.round(-this.m_ScrollPix/ppb);
            from = to - view_seq_width + 1;
        } else {
            from = this.m_FromSeq + Math.round(-this.m_ScrollPix/ppb);
            to = from + view_seq_width -1;
        }

        if (to >= this.m_App.m_SeqLength) to = this.m_App.m_SeqLength - 1;

        return [from, to, to - from + 1];
    },

    createPanel: function() {
        var modeButtHandler = function(b) {
            if (b.iconCls.indexOf('xsv-mode') == 0) {
                this.m_slimMode = b.pressed;
                b.setIconCls(b.pressed ? 'xsv-mode_orig' : 'xsv-mode_slim');
                b.setTooltip('Switch to ' + (b.pressed ? 'normal' : 'slim') + ' mode');
                this.pingClick('2-2-11');
            } else {
                this.setGeneMode(b.pressed);
                this.pingClick('2-2-12');
            }
            this.refresh();
        };
        var slimModeButton = {iconCls: this.m_slimMode ? 'xsv-mode_orig' : 'xsv-mode_slim', scope: this,
            tooltip: 'Switch to ' + (this.m_slimMode ? 'normal' : 'slim') + ' mode',
            pressed: this.m_slimMode, enableToggle: true, handler: modeButtHandler};
        var geneModeButton = {iconCls: 'xsv-gene_mode', itemId: 'geneModeButton', scope: this,
            pressed: false, enableToggle: true, handler: modeButtHandler};
        var seqInfo = this.m_App.m_Config.SeqInfo;

        this.m_DivId = 'graphical_id' + this.m_Idx;
        this.ping({"embedded": this.m_App.m_Embedded});
        if (this.m_App.m_Embedded && this.m_App.m_Embedded == 'minimal') {
            this.m_View = this.m_App.addView({
                layout:'hbox',
                header: false,
                border: false,
                align: 'stretchmax',
                items: [{
                    xtype: 'buttongroup',
                    cls: 'hidden_for_print',
                    columns: 1,
                    width: 24,
                    frame: false,
                    items: [
                        {iconCls: 'xsv-properties', tooltip:'Open Entrez View', scope: this,
                         handler:function(b,e) {this.openFullView(); this.pingClick('2-2-2-1');}},
                        {iconCls:'xsv-flip-strands', id:'flip-button-'+this.m_Idx, tooltip:'Flip Sequence Strands',
                         pressed: this.getFlip(), enableToggle:true, hidden: !this.canFlip(),
                         scope: this, handler:function() { this.flipStrand(); this.pingClick('2-2-7-2');;}},
                        {iconCls: 'xsv-zoom_plus', tooltip:'Zoom In', scope:this,
                         handler:function() {this.zoomIn(); this.pingClick('2-2-5-1');}},
                        {iconCls: 'xsv-zoom_minus', tooltip:'Zoom Out', scope:this,
                         handler:function() {this.zoomOut(); this.pingClick('2-2-5-1');}},
                        {iconCls: 'xsv-zoom_seq', tooltip:'Zoom To Sequence', scope:this,
                         handler:function() {this.zoomSeq(); this.pingClick('2-2-6');}},
                        {iconCls:'xsv-goto_position', tooltip:'Go To Position/Range', scope: this,
                         handler:function() {this.gotoPositionDlg(); this.pingClick('2-2-3');}},
                        geneModeButton,
                        slimModeButton,
                        {iconCls:'x-tbar-loading', tooltip:'Reload View', scope: this,
                         handler:function() {this.refresh(); this.pingClick('2-2-9');}},
                        {iconCls: 'xsv-config', tooltip:'Configure tracks', scope: this, hidden: (this.m_App.m_NoConfDlg === true || this.m_App.m_PermConfId),
                         handler:function() { this.m_App.showTracksConfigDlg('2-2-8-0'); }},
                        {iconCls: 'xsv-question', tooltip:'Help', scope: this,
                         handler:function() {SeqView.showHelpDlg(); this.pingClick('2-2-10');}}
                    ]
                },{
                    header: false,
                    border: false,
                    view:this,
                    flex: 1,
                    autoHeight: true,
                    layout: 'fit',
                    html: {
                        tag: 'div', id: this.m_DivId,
                        cls: 'graphical_div sv-drag sv-highlight sv-dblclick'
                    }
                }]
            });
        } else {
            var tbar = [];
            if (this.m_App.m_Embedded === false) {
                this.m_Spacer = this.m_App.addView({
                    border: false,
                    height: this.getSpacerHeight()

                });
                if (!this.m_Color || this.m_Color.length == 0) {
                    var cp = new Ext.ColorPalette();
                    this.m_Color = cp.colors[Math.floor( Math.random() * cp.colors.length)];
                }
            }

            var range = this.toSeq();
            var r_range_from = this.m_App.posToLocal(range[0]);
            var r_range_to   = this.m_App.posToLocal(range[1]);
            if( this.getFlip() ){
                var t = r_range_from;
                r_range_from = r_range_to;
                r_range_to = t;
            }

            if( this.m_App.m_Toolbar["history"] == true ){
                tbar.push(
                    {iconCls: 'back', tooltip: 'Go back', itemId: 'histPrev', scope:this,
                     handler:function() {this.stepHistory(false, {from_ui: true}); this.pingClick('2-2-1');}}
                );
            }

            if( this.m_App.m_Toolbar["name"] == true ){
                var tiptitle = '<b>' + this.m_App.m_ViewParams['id'] + ': '
                + this.m_App.m_Config.SeqInfo.title + '</b><br><br>' + r_range_from.commify()
                + '&nbsp;-&nbsp;' + r_range_to.commify() + '&nbsp;('
                + (range[1] - range[0] + 1).commify() + '&nbsp;'
                + (this.m_App.m_ViewParams['acc_type']=='protein' ? 'residues' : 'bases');
                + '&nbsp;shown';

                if (this.canFlip())
                    tiptitle += ',&nbsp;' + (this.getFlip() ? 'negative' : 'positive') + 'strand';
                tiptitle += ")";

                var menu_items = [
                    {text: "Entrez View", tooltip: "Switch Entrez View", scope:this, iconCls: 'xsv-full_view',
                     handler: function() {this.openFullView(); this.pingClick('2-2-2-1');}},
                    {text:'Full View', tooltip:'Switch to full standalone view', scope: this,
                     handler: function() {this.openFullView('full'); this.pingClick('2-2-2-2');}}
                ];
                if (seqInfo.assm_info) {
                    menu_items.push({itemId: 'assmInfo', text: 'Assembly ' + seqInfo.assm_info.name, scope: this,
                        tooltip: 'Assembly context information',
                        icon: seqInfo.icon,
                        handler: function() { this.showAssmInfo(); this.pingClick('2-2-2-3'); }
                    });
                }
                

                if (seqInfo.placements && seqInfo.placements.length) menu_items.push('-');
                Ext.each(seqInfo.placements, function() {
                    var link = SeqView.base_url
                        + '?id=' + this.accession
                        + '&v=' + (this.from + 1) + ':' + (this.to + 1);
                    menu_items.push({
                        text: this.label || this.accession + ': ' + (this.from + 1) + '..' + (this.to + 1),
                        tooltip: this.help,
                        handler: function() { window.open(link); }
                    });
                });

                tbar.push({
                    text: 'title',
                    tooltip: tiptitle,
                    itemId: 'tbtitle',
                    menu: new Ext.menu.Menu({ items: menu_items }),
                    scope: this
                });
            }
            if (this.m_App.m_Toolbar["search"] == true) {
                var tt = '<b>Go to position/range:</b><br>'
                        + 'Range formats are 10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m <br><br>'
                        + '<b>Search:</b><br>'
                        + 'Feature name, component name, HGVS, SNP rs id, track info, or sequence (nucleotide regexp with IUPAC equivalents, PROSITE patterns)';
                var buttonSGT = new Ext.Button({
                    text : 'Find:',
                    itemId: 'btnFind',
                    scope: this,
                    handler: function() {this.gotoAndSearch(); this.pingClick('2-2-3');},
                    tooltip: tt
                });
                var gotoBox = Ext.create('Ext.form.field.ComboBox', {
                    itemId: 'gotoBox', width:150,
                    hidden: true, queryModel: 'local', store: this.m_App.searchPatternData,
                    valueField: 'pattern', displayField: 'pattern', typeAhead: true,
                    listeners: {
                        specialkey: function(f,e){
                            if (e.getKey() == e.ENTER){
                               this.gotoAndSearch();
                               this.pingClick('2-2-3');
                               e.stopEvent();//to stop propagation affecting other non-SV text fields in the embedded sviewer
                            }
                        }.bind(this)
                    }
                });

                tbar.push(this.m_App.m_ViewParams.acc_type != 'DNA'? {hidden:true} : '-', buttonSGT, gotoBox);
                //smaller footprint of the above "Find on Sequence" button and its text field for cases when browser's width is
                //less than 900 px
                var buttonSearch = new Ext.Button({
                    itemId: 'btnSearch',
                    hidden: true,
                    tooltip: tt,
                    iconCls: 'xsv-search-button',
                    handler: function() { this.m_App.showSearchParamsDlg(this); this.pingClick('2-2-3');},
                    scope: this
                });
                tbar.push(buttonSearch,'-');
            }

            if( this.m_App.m_Toolbar["panning"] == true ){
                tbar.push(
                    {iconCls:'pan-left', tooltip: 'Pan left', scope:this,
                     handler:function() {this.scrollViewTo( this.m_ScrollPix + 100, SeqView.PAN_LEFT); this.pingClick('2-2-4');}},
                    {iconCls:'pan-right', tooltip: 'Pan right', scope:this,
                     handler:function() {this.scrollViewTo( this.m_ScrollPix - 100, SeqView.PAN_RIGHT); this.pingClick('2-2-4');}},
                    '-'
                );
            }

            if( this.m_App.m_Toolbar["zoom"] == true ){
                tbar.push(
                    {iconCls:'xsv-zoom_minus', tooltip:'Zoom Out', scope: this,
                     handler: function() {this.zoomOut(); this.pingClick('2-2-5-1');}}
                );

                this.m_TbSlider = new Ext.Slider({
                    name: 'zoom', width: 100,
                    minValue: 0, maxValue: 100, value: 0, increment: 5,
                    topThumbZIndex: 100,
                    tipText: function(thumb){ return String(thumb.value) + '%'; },
                    listeners: {
                        changecomplete: function(el,val) {
                            this.pingClick('2-2-5');
                            var slider_range = 100;
                            var slider_val = 100 - val;

                            var vis_range = this.toSeq();
                            var zoom_point = vis_range[0] + vis_range[2] / 2;

                            var max_bases_to_show = this.m_App.m_SeqLength;
                            var min_bases_to_show = this.getMostRightPix() * SeqView.MinBpp;

                            var ee = Math.log(max_bases_to_show / min_bases_to_show);

                            var slider_val = Math.min(slider_range, Math.max(0, slider_val));

                            var len = min_bases_to_show * Math.exp(ee * slider_val / slider_range);
                            var from = zoom_point - len / 2;

                            len = Math.min(this.m_App.m_SeqLength, Math.round(len));
                            from = Math.round( Math.max(0, from) );
                            if (from + len > this.m_App.m_SeqLength) {
                                var extra = from + len - this.m_App.m_SeqLength;
                                from -= extra;
                            }
                            this.startImageLoading(from, len, {from_ui: true});

                        }.bind(this)
                    }
                });
                tbar.push( this.m_TbSlider );

                tbar.push(
                    {iconCls:'xsv-zoom_plus', tooltip: 'Zoom In', scope: this,
                     handler:function() {this.zoomIn(); this.pingClick('2-2-5-1');}},
                    {iconCls: 'xsv-zoom_seq', tooltip: 'Zoom To Sequence', scope: this,
                     handler:function() {this.zoomSeq(); this.pingClick('2-2-6');}}
                );
            }

            if (this.m_App.m_Toolbar["genemode"] == true) tbar.push(geneModeButton);

            tbar.push('->',
                {id:'seq-view-sync' + this.m_Idx, hidden: this.m_App.m_NeedAlignmentView == false, text: 'Sync Alignment View', scope: this,
                 handler: function(){this.m_App.onSyncAlignView(this);} },
                {iconCls: 'xsv-full_view', tooltip:'Open Entrez View', scope: this, hidden: (this.m_App.m_Embedded === true || this.m_App.m_NoViewHeader !== true),
                 handler: function(b,e){this.openFullView();}}
            );

            if (this.m_App.m_Toolbar["tools"] == true) {
                var menu_items = [
                    {text:'Go To', iconCls:'xsv-goto_position', tooltip:'Go To Position/Range', scope: this,
                     handler:function(){this.gotoPositionDlg(); this.pingClick('2-2-7-0');}},
                    {text:'Search', iconCls:'xsv-search-button', tooltip:'Search features, components or sequences', scope: this,
                     handler:function(){this.m_App.showSearchParamsDlg(this); this.pingClick('2-2-7-1');}},
                    {text:'Flip Strands', tooltip:'Flip Sequence Strands', iconCls:'xsv-flip-strands',
                     pressed: this.getFlip(), enableToggle:true, hidden: !this.canFlip(), scope: this,
                     handler:function(){this.flipStrand(); this.pingClick('2-2-7-2');}},
                    {text:'Markers', iconCls:'xsv-markers', scope:this,
                     handler:function(){this.m_App.showMarkersDlg(this); this.pingClick('2-2-7-3');}},
                    {text:'Set Origin', iconCls:'xsv-origin', scope:this,
                     handler:function(){this.m_App.showOriginDlg(null); this.pingClick('2-2-7-4');}},,
                    {text:'Sequence Text View', iconCls:'xsv-new_fasta', tooltip:'Create New Sequence Text View', scope: this,
                     handler:function(){this.m_App.createTextView(this); this.pingClick('2-2-7-5');}}
                ];

                if( this.m_App.mf_MultiPanel ){
                    menu_items.push(
                        '-',
                        {text:'Add new Panel', iconCls:'xsv-new_view', tooltip:'Create New Graphical Panel', scope:this,
                        handler:function() {
                            var view = new SeqView.Graphic(this.m_App);
                            this.m_App.registerView(view);
                            this.pingClick('2-2-7-6');}
                        },

                        {text:'Add new Panel on Range', iconCls:'xsv-new_view', tooltip:'Create New Graphical Panel on Selected Range',
                            scope:this, handler:function() {
                                var view = new SeqView.Graphic(this.m_App);
                                this.m_App.registerView(view);
                                var range = this.getTotalSelectedRange();
                                if (range) this.removeRangeSelection();
                                else
                                    range = this.m_UrlFrom ? [this.m_UrlFrom, this.m_UrlTo] :
                                        [this.m_VisFromSeq + 1, this.m_VisFromSeq + this.m_VisLenSeq];

                                view.startImageLoading(range[0], range[1] - range[0] + 1, {from_ui: true} );
                                this.pingClick('2-2-7-7');
                            }
                        },
                        {text: 'Show/Hide Status bar', scope: this, handler: function() {
                            var sb = this.m_View.down('#statusBar');
                            sb[sb.isVisible() ? 'hide' : 'show']();
//                            localStorage.setItem('NCBI/SV/statusBar', sb.isVisible() ? 'hide' : 'show');
//                            sb[localStorage.getItem('NCBI/SV/statusBar')]();
                        }}
                    );
                }
                var noPBlast = (this.m_App.m_ViewParams['acc_type'] !== 'DNA');
                menu_items.push(
                    '-',
                    {text:'BLAST and Primer Search', menu:new Ext.menu.Menu({items:[

                    {text:'BLAST Search (Visible Range)', iconCls:'xsv-blast', scope: this,
                     handler: function() {
                        this.m_App.blast([this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1]);
                        this.pingClick('2-2-7-8-0');}},
                    {text:'Primer BLAST (Visible Range)', iconCls:'xsv-primer', scope: this, disabled: noPBlast,
                     handler: function() {
                        this.primerBlast(false, [ [this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1] ]);
                        this.pingClick('2-2-7-8-1');}},
                    {text:'BLAST Search (Selection)', iconCls:'xsv-blast', scope: this,
                     handler: function() {this.blastSelection(); this.pingClick('2-2-7-8-2');}},
                    {text:'Primer BLAST (Selection)', iconCls:'xsv-primer', scope: this, disabled: noPBlast,
                     handler: function() {this.primerBlast(); this.pingClick('2-2-7-8-3');}}
                    ]})},
                    '-',
                    {text:'Download', iconCls:'xsv-download-static',
                     menu:new Ext.menu.Menu({items:[
                        {text:'FASTA (Visible Range)', scope: this,
                        handler: function () { this.pingClick('2-2-7-9-0'); this.downloadData(false, "fasta", null); }},
                        {text: 'FASTA (All Markers)', scope: this,
                        handler: function () { this.pingClick('2-2-7-9-5'); this.downloadAllMarkers("fasta"); }},
                        {text:'FASTA (All Selections)', scope: this,
                         handler:function() {this.pingClick('2-2-7-9-1'); this.downloadData(true, "fasta", null);}},
                        {text: 'GenBank Flat File (Visible Range)', scope: this,
                         handler: function() {this.pingClick('2-2-7-9-2'); this.downloadData(false, "flat", null);}},
                        {text:'GenBank Flat File (All Selections)', scope: this,
                         handler:function() {this.pingClick('2-2-7-9-3'); this.downloadData(true, "flat", null);}},
                        {text:'PDF file', disabled: this.m_App.m_NoPDF, scope: this,
                         handler:function() {this.pingClick('2-2-7-9-4'); this.downloadPDF();}}
                    ]})},
                    '-',
                    {text:'Printer-Friendly PDF', iconCls:'xsv-printer', disabled: this.m_App.m_NoPDF, scope: this,
                        handler: function() {this.pingClick('2-2-7-10'); this.downloadPDF();}},
                    {text: 'Preferences', iconCls: 'xsv-configure', scope: this,
                        handler: function() {this.pingClick('2-2-7-11'); this.setPreferences();}}
                );

                var tools_menu = new Ext.menu.Menu({ items: menu_items });

                tbar.push({
                    text: 'Tools', iconCls: 'xsv-tools', tooltip: 'Tools', itemId: 'buttonTools', scope: this,
                    handler: function() {this.pingClick('2-2-7');},
                    menu: tools_menu
                 });
            }

            if (this.m_App.m_Toolbar["slimmode"] == true) tbar.push(slimModeButton);

            var cfgButt = {text: 'Tracks', iconCls: 'xsv-config', tooltip:'Configure tracks', scope: this,
                           handler: function(e) { this.m_App.showTracksConfigDlg('2-2-8-0');  }};
            var cfgPanel = (this.m_App.m_NoConfDlg !== true
                         && !this.m_App.m_PermConfId
                         && this.m_App.m_Toolbar["config"] == true);
            this.m_App.m_Toolbar["trackSets"] = (TMS.TrackSets != undefined && this.m_App.m_AppContext != null);
            if (this.m_App.m_Toolbar["trackSets"] == true) {
                var icon = 'xsv-search-loading';
                var stdID = 'sv_stdtrackset_id' + this.m_Idx + '_' + this.m_App.m_Idx;
                var menu_items = [];
                if (cfgPanel) {
                    menu_items.push(
                        {text: 'Configure tracks', tooltip: 'Configure current set of tracks', scope: this,
                        iconCls: 'xsv-config', handler: cfgButt.handler}, '-');
                }
                cfgButt.handler = function(butt) {
                    if (butt.menu.hidden) return;
                    this.pingClick('2-2-8');
                    if (!SeqView.requestTrackSets) {
                        delete this.m_App.m_currentTrackSetId;
                        SeqView.configureTrackSets(this.m_App);
                    }
                    var cmp = Ext.getCmp(stdID).disable();
                    cmp.setIconCls(icon);
                    cmp.menu.removeAll();
                    cmp = cmp.nextSibling().disable();
                    cmp.setIconCls(icon);
                    cmp.menu.removeAll();
                    var view = this;
                    SeqView.requestTrackSets(function(s) { view.processTrackSets(s, stdID); });
                };
                menu_items.push({text: 'NCBI Recommended Track Sets', id: stdID,
                                iconCls: icon, menu: new Ext.menu.Menu({items: []}), disabled: true},
                                {text: 'My NCBI Track Collections',
                                iconCls: icon, menu: new Ext.menu.Menu({items: []}), disabled: true},
                                {text: 'FAQ', scope: this,
                                iconCls: 'xsv-question', handler: function() {
                                    this.pingClick('2-2-8-3', 'goto_FAQ_page');
                                    window.open(SeqView.webNCBI + 'tools/sviewer/faq/#tracksets');}});
                cfgButt.menuAlign = 'tr-br?';
                cfgButt.menu = new Ext.menu.Menu({ items: menu_items });
                tbar.push('-', cfgButt);
            } else if (cfgPanel) tbar.push('-', cfgButt);

            if (this.m_App.m_Toolbar["reload"] == true){
                tbar.push(
                    {iconCls:'x-tbar-loading', disabledCls: 'xsv-search-loading',
                     id: 'seq-view-loading-' + this.m_Idx, tooltip: 'Reload View',
                     scope: this, handler:function(){this.refresh(); this.pingClick('2-2-9');}}
                );
            }

            var first_graphic = true;
            for (var i = 0; i < this.m_App.m_Views.length; i++) {
                var view = this.m_App.m_Views[i];
                if (view && view != this && view.isGraphic()) {
                    first_graphic = false;
                    break;
                }
            }
            if (first_graphic) {
                if( this.m_App.m_Toolbar["help"] == true ){
                    tbar.push(
                        {iconCls: 'xsv-question',
                         tooltip:'Help',
                         layout: {type:'vbox'},
                         scope: this,
                         handler: function() {this.pingClick('2-2-10');},
                         menu: new Ext.menu.Menu(
                            {defaultOffsets: [-70,0], scope: this,
                             items: [
                                {text: 'Help', iconCls: 'xsv-question', scope:this,
                                 handler: function() {SeqView.showHelpDlg(); this.pingClick('2-2-10-1');}},
                                {text: 'Link to View', iconCls: 'xsv-link_to_page', scope: this,
                                 handler: function() {this.m_App.showLinkURLDlg(); this.pingClick('2-2-10-3');}},
                                {text: 'Feedback', iconCls: 'xsv-feedback', scope: this,
                                 handler: function() {this.m_App.showFeedbackDlg('2-2-10-0');}},
                                {text: 'About', scope: this,
                                 handler: function() {SeqView.showAboutMessage(this.m_Idx); this.pingClick('2-2-10-2');}}
                            ]})
                        }
                    );
                }
            } else {
                tbar.push(
                    {iconCls: 'tb-close', tooltip:'Close', scope: this, handler:function() { this.remove(); }}
                );
            }

            var tools = [];
            if (this.m_App.m_Embedded === false) {
                tools.push(this.createCloseTBar());
            } else {
                tools.push({id:'help',qtip: 'Help', handler: function(){this.pingClick('2-2-10-1'); SeqView.showHelpDlg();}})
            }
            this.m_View = this.m_App.addView({
                collapsible: true, title: 'Loading...', header: false,
                tools:tools, view:this, tbar:tbar,
                dockedItems: [{ itemId: 'statusBar', 
                    xtype: 'toolbar', dock: 'bottom', //hidden: localStorage.getItem('NCBI/SV/statusBar') != 'show',
                    items: [
                        {itemId: 'status', xtype: 'tbtext', text: ''},
                        '->',
                        {itemId: 'feedback', xtype: 'button', text: '', scope: this.m_App, iconCls: 'xsv-feedback',
                        tooltip: 'Feedback', 
                        handler: function() { this.m_App.showFeedbackDlg('10-2'); }},
                        {itemId: 'tracks', xtype: 'button', text: '', scope: this.m_App, iconCls: 'xsv-config',
                        handler: function () { this.showTracksConfigDlg('10-1'); }}

//                        {itemId: 'link', xtype: 'tbtext', text: ''}
                    ]
                }],
                html: {
                    tag: 'div', id: this.m_DivId,
                    cls: 'graphical_div sv-drag sv-highlight sv-dblclick'
                }
            })
        }
        var the_div = Ext.get(this.m_DivId);
        this.m_vLine = the_div.appendChild(new Ext.Element(document.createElement('div')));
        this.m_vLine.set({id: 'lRuler' + this.m_Idx,
             style: 'position:absolute;border:1px solid lightgray;border-right-width:0;top:0;width:0;z-index:10;display:none;'});
/*if (Ext.supports.Touch) {
        the_div.on({ scope: this,
            pointerdown:  this.onMouseDown,
            pointerup:   this.onMouseUp,
            dblclick:   this.onDblClick,
            pointermove:  this.onMouseMove,
            contextmenu:this.onContextMenu,
            pointerout:   function() {this.m_vLine.setStyle('display', 'none');}
        });
} else */{
        the_div.on({ scope: this,
            mousedown:  this.onMouseDown,
//            tap:        this.onTap,
            touchstart: this.onMouseDown,
            touchend:   this.onMouseUp,
            dblclick:   this.onDblClick,
            mousemove:  this.onMouseMove,
            touchmove:  this.onMouseMove,
            contextmenu:this.onContextMenu,
            doubletap:  this.onContextMenu,
            mouseout:   function() {this.m_vLine.setStyle('display', 'none');}
        });
}
        if (seqInfo.warnTT) {
            var btn = this.m_View.down('#assmInfo');
            if (btn) btn.setTooltip('WARNING!<br>' + seqInfo.warnTT);
        }
        this.updateStatusBar('feedback', '');
        this.startImageLoading(this.m_FromSeq, this.m_LenSeq, {fire_event:false});
//        this.updateStatusBar('link', '<a href=\"https://www.ncbi.nlm.nih.gov/tools/sviewer/\" target=\"_blank\">Help</a>');
    },
    onTap: function(e) {
        //if (e.type == 'touchstart')
        {
            console.log(e.type);
            e.stopEvent();
        }
    },

    updateStatusBar: function(prop, text, icon) {
        var st = this.m_View.down('#' + prop);
        if (st) {
            if (typeof text == 'string') st.setText(text);
            if (typeof icon == 'string') st.setStyle('background', icon + (icon ? ' no-repeat' : ''));
            st.setStyle('padding', (st.el.getStyle('background-image') != 'none' ? '0 22px' : '0 3px'));
        }
        return st || {hide: Ext.emptyFn, show: Ext.emptyFn};
    },

    processTrackSets: function(sets, stdID) {
        var scope = this;
        var prepareTT = function() {
            var tt = this.tt;
            if (!tt) {
                if (this.type == 'default') {
                    this.tracks = scope.m_App.m_defaultTrackSet || [];
                }
                Ext.each(this.tracks, function(t, i) {
                    tt += '<p style="margin:0;' + (i%2 ? '">' : 'background-color:gainsboro;">') + t.GetName() + '</p>';
                });
            }
            if (tt) this.tooltip = new Ext.ToolTip({target: this.el, title: this.text, html: tt,
                anchor: 'left', margin: 0, anchorOffset: 185, track_Mouse: true, dismissDelay: 0});
        }
        var usrSets = [];
        var tmsSets = [{text: 'Reset to Default Tracks',
                        type: 'default', tsId: 'defaultTS', tt: '',
                        handler: function() { scope.applyTrackSet(this); },
                        listeners: { afterrender: prepareTT }}];
        Ext.each(sets, function() {
            var tSet = {text: this.GetName(),
                        tracks: this.GetTracks(), tt: '', tsId: this.GetId(),
                        handler: function() { scope.applyTrackSet(this); },
                        listeners: { afterrender: prepareTT }};
            switch (this.GetTrackSetType()) {
                case TMS.TrackSetType.TMS: tmsSets.push(tSet); tSet.type = 'tms'; break;
                case TMS.TrackSetType.myNCBI_Collection: usrSets.push(tSet); tSet.type = 'MyNCBI'; break;
            }
        });

        var cmp = Ext.getCmp(stdID).enable();;
        var addMenuItem = function() {
            if (scope.m_App.m_currentTrackSetId == this.tsId) this.iconCls = 'xsv-search-results';
            cmp.menu.add(this);
        }
        cmp.setIconCls('xsv-seq-logo');
        Ext.each(tmsSets, addMenuItem);
        cmp = cmp.nextSibling().enable();
        cmp.setIconCls('xsv-seq-logo');
        if (usrSets.length) {
            Ext.each(usrSets, addMenuItem);
            cmp.menu.add('-');
        }
        var inMyNCBI = MyNCBI.User.IsLoggedIn();
        cmp.menu.add({text: !inMyNCBI ? 'Sign in to NCBI' : 'Manage collections', scope: this,
                handler: function(){
                    this.pingClick('2-2-8-2', inMyNCBI ? 'NCBI_sign_in' : 'collections_manage');
                    window.open(SeqView.webNCBI + 'myncbi/collections');}});
        if (inMyNCBI)
            cmp.menu.add({text: 'Save current tracks', disabled: this.m_App.m_noSaveTracks, handler: function() { scope.saveTrackSet(); }});
    },

    saveTrackSet: function() {
        var idList = [],
            dOpts = {};

        this.pingClick('2-2-8-1', 'saveTrackSet');
        Ext.each(this.m_App.getActiveTracks(), function() {
            dOpts[this.id] = SeqView.TM.getTrackDisplayOptions(this).slice(0);
            idList.push(this.id);
        });
        if (!idList.length) return;
        TMS.GetTracksById(idList).done(function(){
            var tracklist = this;
            Ext.MessageBox.prompt('Save Track Collection', 'Please enter collection name:', function(btn, name) {
                if (btn!='ok'  || name.length==0) return;
                var tracks = [];
                Ext.each(tracklist.GetTracks(), function() {
                    var dti = new TMS.DisplayTrackInfo(this);
                    dti.SetDisplayOptions(dOpts[this.GetTMSId()]);
                    tracks.push(dti);
                });

                TMS.TrackSets.TrackSetService.CreateTrackset(name, tracks).done(function(trackset_id) {
                    console.log('trackset ' + trackset_id + 'has been added');
                });

            });
       });
    },

    applyTrackSet: function(tSet) {
        this.ping({'jsevent': 'click', 'sv-area': '2-2-8-' + tSet.type,
                   'sv-event': 'applyTrackSet', trackSetID: tSet.tsId});
        var trx = '',
            ts = {tracks: []};
        Ext.each(tSet.tracks, function() {
            ts.tracks.push({displayOptions: this.GetDisplayOptions(), id: this.GetTMSId()});
        });
        this.m_App.fireEvent('trackSet_requested', ts);
        Ext.each(ts.tracks, function() {
            trx += '[id:' + this.id + (this.displayOptions ? ',' + this.displayOptions : '') + ']';
        });
        this.loadTrackSet(trx || '[amend]', tSet.tsId, ts.ignoreURL);
    },

    loadTrackSet: function(tracks, tsId, ignoreURL) {
        var app = this.m_App;
        this.setGeneMode(false);
        var dflt = app.defaultConfig.m_TracksFromURL || '';
        if (ignoreURL) dflt = '';

        if (tsId !== 'defaultTS')
            app.m_TracksFromURL = dflt.replace(/]/g, ',shown:false]') + tracks;
        else
            app.m_TracksFromURL = tracks + dflt;
        app.m_currentTrackSetId = tsId;
        var load_cmp = Ext.getCmp('seq-view-loading-' + this.m_Idx).disable();

        app.m_Config.load({forcereload: true,
            callback: { scope: this,
                success: function(res){ app.fireEvent('configuration_changed', app); this.refresh(); },
                failure: function(data, text, res) { load_cmp.enable(); app.infoMissed(data, text, res); }
            }
        });
    },

    clear: function() {
        this.m_History = [];
        this.m_HistIx = 0;
        this.mf_HistoricalBack = false;

        this.m_Selection = null;
        this.m_RangeSelectionSet = [];
        this.m_Reflections = null;
        this.m_SelectedSig = null;

        if( this.m_App){
            this.m_App.fireEvent( 'selection_changed', this );
        }

        this.m_ScrollPixBeforeLoad = 0;
        this.m_Loading = false;

        SeqView.Graphic.superclass.clear.call(this);
    },

//////////////////////////////////////////////////////////////////////////
// refresh:

    refresh: function(options) {
        if (this.m_VisLenSeq == 0)
            return;
        this.startImageLoading(this.m_VisFromSeq, this.m_VisLenSeq, options);
    },

    // remove hover selection decoration
    removeFloatingSelection: function() {
       if (this.m_Selection) {
           this.m_Selection.remove();
           this.m_Selection = null;
       }
    },

    selectRange: function(sel_range) { this.selectRangeSet([sel_range]); },

    selectRangeSet: function(range_set, noEvent) {
        if (!range_set || !range_set.length) return;
        this.removeRangeSelection(null, noEvent);
        for (var x in range_set) {
            this.m_RangeSelectionSet.push(new SeqView.RangeSelection(this, range_set[x]));
        }
    },
    //removes RangeSelections and SelectedRangeSets when no pinned range tooltip associated with them
    removeSelectionsWithNoPinnedToolTips: function() {
        var length = this.m_RangeSelectionSet.length;
        for (var x = length - 1; x >= 0; x--) {
            var RSSx = this.m_RangeSelectionSet[x];
            //!!!RSSx.m_ToolTip.remove(); // ???
            RSSx.remove();
            this.m_RangeSelectionSet.splice(x, 1);
        }
        if (length != this.m_RangeSelectionSet.length) this.m_App.fireEvent('selection_changed', this);
    },

    popRangeElement: function(range_elt) {
        var sel_set = this.m_RangeSelectionSet;
        for (var x in sel_set) {
            if (sel_set[x] == range_elt) {
                sel_set.splice(x, 1);
                break;
            }
        }
    },

    removeRangeSelection: function(range, noEvent) {
        var sel_set = this.m_RangeSelectionSet;
        if (range) {
            for (var x in sel_set) {
                if (sel_set[x].range != range) continue;
                sel_set.splice(x, 1)[0].remove();
                break;
            }
        } else {
            Ext.each(sel_set, function(rs) { rs.remove(); });
            this.m_RangeSelectionSet = [];
        }
        if (!noEvent) this.m_App.fireEvent('selection_changed', this);
    },

    removeFeatMarks: function() {
        if (this.m_HighLightedFeature) {
            Ext.each(this.m_HighLightedFeature,function(sf) {
                sf.remove();
            });
        }
        delete this.m_HighLightedFeature;
    },

    startRangeSelection: function(pix_start, ctrl) {
        // this object is for user interaction time only, after
        // selection is finished, it will be properly registered
        // or discarded, see addRangeSelection
        this.m_InDragAction = true;
        new SeqView.RangeSelection(this, [pix_start], true, ctrl);
    },

    addRangeSelection: function(sel_element, range, pixel_width) {
        this.m_InDragAction = false;
        if (pixel_width < 5) {
            sel_element.remove();
            return false;
        }
        sel_element.updateCoords();
        this.m_RangeSelectionSet.push(sel_element);
        this.m_App.fireEvent( 'user_changed_selection', this );
        return true;
    },

    getTotalSelectedRange: function() {
        if (!this.m_RangeSelectionSet.length) return null;
        var min_pos = this.m_App.m_SeqLength;
        var max_pos = 0;
        Ext.each(this.m_RangeSelectionSet, function() {
            var rng = this.range;
            max_pos = Math.max(max_pos, rng[0], rng[1]);
            min_pos = Math.min(min_pos, rng[0], rng[1]);
        });
        return [min_pos, max_pos];
    },

    getSelectedRangeSet: function() {
        var srs = [];
        Ext.each(this.m_RangeSelectionSet, function(){ srs.push(this.range); });
        return srs;
    },

    startImageLoading: function(from, len, options) {
        if (!arguments.length) {
            from = this.m_deferLoadParams[0];
            len = this.m_deferLoadParams[1];
            options = this.m_deferLoadParams[2];
        }
        if (this.m_Loading) {
            var info = { target: this, method: this.startImageLoading, args: arguments };
            this.m_App.fireEvent('defer_image_load', info);
            if(!info.ignore) {
                if (!this.m_deferLoadParams) 
                   Ext.defer(this.startImageLoading, 1000, this);
                if (arguments.length) this.m_deferLoadParams = arguments;
            }
            return;
        };
        delete this.m_deferLoadParams;
        this.m_Loading = true;
        this.clearCache();
        this.loadImage(from, len, options);
    },

    // Check that the range fits sequence boundaries,
    // if not, adjust it appropriately
    adjustRequest: function(from, len) {
        if (from < 0) from = 0;
        if (from + len > this.m_App.m_SeqLength) {
            var over = from + len - this.m_App.m_SeqLength;
            from -= over;
            if (from < 0) {
                from = 0;
                len = this.m_App.m_SeqLength;
            }
        }
        return [from, len];
    },


    preinitViewParams: function(data) {
        if (this.m_bpPix) return;
        this.m_Width = data.img_width;
        this.m_LenSeq = data.len;
        this.m_FromSeq = data.from;
        this.m_bpPix = (this.m_Width || SeqView.ChunkWidth)/this.m_LenSeq;
        
        var screen_width = this.getScreenWidth(); //the_div.getWidth();
        var off_pix = Math.min(Math.max(0,this.m_Width - screen_width), this.seq2Pix(this.m_ReqFrom-0.5));
        if (this.getFlip()) {
            off_pix = Math.max(0, this.seq2Pix(this.m_ReqFrom - 0.5) - screen_width);
        }

        // adding scrolling happend during request
        this.m_ScrollPix -= this.m_ScrollPixBeforeLoad ? this.m_ScrollPixBeforeLoad : this.m_ScrollPix; 
        this.m_ScrollPixBeforeLoad = 0;

        this.m_ScrollPix -= off_pix;
    },


    requestImages: function(url, params, reqNum, options) {
        this.m_App.m_timeStamps.loadImagesTS = new Date().getTime();
        this.m_Height = this.m_bpPix = 0;
        this.tmp_showGeneModeButton = false;
        this.graphTracksNew = [];
        this.m_FromCgi = {areas:[], images:[]};
        this.updateStatusBar('status', 'Loading...', SeqView.ExtIconLoading);
        this.updateStatusBar('tracks').hide(); 

        var app = this.m_App,
            gview = this,
            tmp_its = [],
            imgRQs = params.tracks,
            numRQs = params.tracks.length;

        var processResponse = function(data, idx) {
            if (reqNum < this.m_ReqNum) return;
            if (!data.job_status) {
                var it = gview.m_FromCgi.images[idx] = {img: new Image(), idx: idx, h: data.img_height, data_height: data.data_height, areas: data.areas};
                gview.preinitViewParams(data);
                gview.m_FromCgi.areas = gview.m_FromCgi.areas.concat(data.areas || []);

                if (data.hairlines) {
                    it.hairlines = data.hairlines;
                    data.hairlines.forEach(function(hl){hl.pos = Math.floor(hl.pos);});
                }
                var img_url = data.img_url;
                if (data.img_url && data.img_url.charAt(0) == '?') {
                    img_url = app.m_CGIs.NetCache + data.img_url;
                }

                it.img_url = it.img.src = img_url;

                it.img.onload = function() {
                    tmp_its[idx] = it;
                    if (SeqView.scrolled && numRQs > 1) return numRQs--;
                    for (var x in tmp_its) { if (!NCBIGBUtils.isNumeric(x)) continue;
                        gview.addTracksImage(tmp_its[x]);
                        delete tmp_its[x];
                    }
                    if (--numRQs) return;
                    gview.processTracksImages();
                    gview.applyResultsOnLoad(options);

                    var warn = (gview.m_App.m_appWarning || '') + gview.m_App.m_Config.SeqInfo.warnTT;
                    delete gview.m_App.m_appWarning;

                    var st = gview.updateStatusBar('status', warn ? 'Warning' : '',
                        warn ? ('url(' + SeqView.ExtIconLoading.slice(5, -18) + 'shared/warning.gif)') : '');
                    if (warn) {
                        if (!st.tooltip) st.tooltip = Ext.create('Ext.tip.ToolTip', {target: st.id, title: 'Warning'});
                        st.tooltip.setHtml(warn);
                    }
                    else
                        if (st.tooltip) st.tooltip = st.tooltip.destroy(); 
                    gview.updateStatusBar('tracks', 'Tracks shown: ' + gview.m_graphTracks.length + '/' + gview.m_App.m_Config.TrackConfig.length).show();
                };
                it.img.empty = (typeof it.areas == 'undefined');
            } else {
                var jst = data.job_status;
                if (jst == 'submitted' || jst == 'running' || jst == 'pending') {
                    var cfg = {
                        url: app.m_CGIs.Graphic + '?job_key=' + data.job_id,
                        success: function(d, t, r) { processResponse.call(gview, d, idx); },
                        error: function(d, t, r) { if (reqNum == gview.m_ReqNum) gview.loadFailure(t, r); }
                    };
                    Ext.defer(SeqView.App.simpleAjaxRequest, 2000, gview, [cfg]);
                }
            }
        };

        imgRQs.forEach(function(t, idx){
            params.tracks = t;
            app.AjaxRequest({
                url: url, data: params,
                success: function(d, t, r) { processResponse.call(gview, d, idx); },
                error: function(d, t, r) { if (reqNum == gview.m_ReqNum) gview.loadFailure(t, r); }});
        });
    },


    loadImage: function(from, len, options){
        if (this.m_featMarkers){
            this.m_featMarkers.div.parentElement.removeChild(this.m_featMarkers.div);
            delete this.m_featMarkers;
        };

        from = Math.floor(Math.max(0, from));
        len = Math.floor(Math.min(len,(this.m_App.m_SeqLength - from)));

        var screen_width = this.getScreenWidth(); // our panel width in pixels
        if (screen_width <= 0 || !this.m_View.isVisible()) {
            // the view is hidden, just save the visible range.
            // the client will have to call refresh method to load the image
            this.m_VisFromSeq = from;
            this.m_VisLenSeq = len;
            this.m_Loading = false;
            this.m_App.m_InitialLoading = false;
            return;
        }

        var req_from = from;
        var req_len = len;

        // Basepairs can not be wider than 1/SeqView.MinBpp (24 pixels per basepair)
        if( req_len < screen_width*SeqView.MinBpp ){
            var new_req_len = Math.floor( screen_width*SeqView.MinBpp );
            if( new_req_len > this.m_App.m_SeqLength ){
                new_req_len = this.m_App.m_SeqLength;
            }
            // we check for range validity and adjust it later
            req_from = Math.floor( req_from - (new_req_len - req_len)/2 );
            req_len = new_req_len;
        }

        // Check that view is detailed, so we can't get away with fractional
        // pixels per basepair. If there are more than 4 pixels per base, bases
        // are becoming individually visible, so fractional ppb will lead to
        // uneven spacing between base letters
        var ppb = screen_width/req_len;
        if (ppb > 4) {
            ppb = Math.floor(ppb);
            var new_req_len = Math.floor(screen_width/ppb);
            if (new_req_len > this.m_App.m_SeqLength) {
                new_req_len = this.m_App.m_SeqLength;
                ppb = screen_width/new_req_len;
            }
            req_len = new_req_len;
        }
        var from_len = this.adjustRequest(req_from, req_len);
        req_from = from_len[0]; req_len = from_len[1];

        var chunk_width;
        if( req_len > 2000000 ){
            chunk_width = 3*screen_width;
        } else if( req_len > 1000000 ){
            chunk_width = 5*screen_width;
        } else {
            chunk_width = 8*screen_width;
        }
        chunk_width = Math.min(SeqView.ChunkWidth, chunk_width);
        // Actually chunk can be little bigger than maximum.

        var fit_len = Math.floor( chunk_width / ppb );

        from = Math.floor( req_from - (fit_len - req_len)/2 );
        len = fit_len;

        var from_len = this.adjustRequest(from, len);
        from = from_len[0]; len = from_len[1];
        chunk_width = len * ppb;

        if (!isNaN(chunk_width)) {
            chunk_width = Math.floor(chunk_width);
        }
         // added some sanity check before calling seqgraphic.cgi
        if (isNaN(from) || isNaN(len) || isNaN(chunk_width) || chunk_width > 16384) {
            // TODO we need to understand why we got in here.
            // Nothing left than just cancel loading and return.
            this.m_Loading = false;
            return;
        }

        this.m_ReqFrom = req_from;
        this.m_ReqLen  = req_len;
        if (this.m_VisFromSeq != req_from || this.m_VisLenSeq != req_len)
            this.updateTitle(options);

        len = Math.max(len, Math.min(this.m_App.m_SeqLength, 100));
        // go for it!!
        var url = this.m_App.m_CGIs.Graphic;
        var params = this.getGraphicParams(from, len, chunk_width);

        if (options && options.data_changed === true) {
            params.data_changed = 'true';
        }
        
        var load_cmp = Ext.getCmp('seq-view-loading-' + this.m_Idx);
        if (load_cmp) load_cmp.disable();
        document.getElementById(this.m_App.m_DivId).style.cursor = 'progress';

        // remove selection decorations
        this.removeFloatingSelection();
        this.m_SelectedRangeSet = this.m_SelectedRangeSet || this.getSelectedRangeSet();
        this.removeRangeSelection();
        this.removeFeatMarks();
        this.m_App.watchdogStart(url, '', params);

        this.requestImages(url, params, ++this.m_ReqNum, options);

        if (this.m_TbSlider) {
            var slider_range = 100;
            var vis_len = req_len || this.toSeq()[2];
            var min_bases_to_show = this.getMostRightPix() * SeqView.MinBpp;

            var ee = Math.log(this.m_App.m_SeqLength / min_bases_to_show);
            var slider_val = slider_range * Math.log(vis_len / min_bases_to_show) / ee;

            this.m_TbSlider.setValue(100 - slider_val);
        }
    },

    getGraphicParams: function(from, len, width, forPDF) {
        var params = {id: this.m_App.GI, client:'seqviewer', width: width,
                      view_width: this.getScreenWidth(), from:from, len:len};
        Ext.apply(params, this.m_App.m_GraphicExtraParams);

        params.forPDF = forPDF;

        if (this.m_SelectedSig) params.select = this.m_SelectedSig;
        if (this.getFlip()) params.flip = 'true';

        this.addURLParams(params);
        params.markers = this.m_App.m_MarkersInfo.getMarkersData(true, forPDF);
        return params;
    },

    // drawing feature markers
    drawFeatMarkers: function() {
        var div = this.m_featMarkers.div = document.createElement('div');
        var trx = this.m_graphTracks,
            top = this.m_topRulerHeight,
            height = this.m_Height - top - this.m_bottomRulerHeight,
            left = this.m_featMarkers.left;
        div.setAttribute('class', 'sv-drag sv-dblclick');
        div.setAttribute('style', 'position: relative; top:' + top
            + 'px; left:' + (left + this.m_ScrollPix)
            + 'px; width:' + this.m_featMarkers.width
            + 'px; height:' + height + 'px;');
        var tmpl = '<div class="sv-drag sv-dblclick"' + //id="hlbox_' + this.m_DivId  + '"
        ' style="border-width: 1px; border-left-style: solid; position: absolute;';

        this.m_featMarkers.hls.forEach(function(hls){
            var arrHLs = hls.hairlines;
            for (var idx = trx.length - 1; idx > 0; idx--) {
                if (trx[idx].divT == hls.top) break;
            }
            var tHeight = trx[idx].divT - top,
                bTop = tHeight + hls.height,
                bHeight = height - bTop;

            for(var i = 0; i < arrHLs.length; i += 2) {
                 var lHL = arrHLs[i],
                     rHL = arrHLs[i + 1] || lHL,
                     lPos = lHL.pos - left,
                     rPos = rHL.pos - lHL.pos+1;

                if (i == arrHLs.length - 1)
                    rPos = -1;

                lHL.top_color = (lHL.top_color.r == 255) ? 'red' : 'grey';
                if (typeof rHL.top_color == 'object') rHL.top_color = (rHL.top_color.r == 255) ? 'red' : 'grey';
                if (tHeight)
                    div.innerHTML += tmpl + 'width:' + rPos + 'px; border-right-style:' + (rPos >= 0 ? 'solid' : 'none')
                        + ';border-left-color:' + lHL.top_color + '; border-right-color:' + rHL.top_color
                        + ';opacity: 0.4;top: 0px; height:' + tHeight + 'px; left:' + lPos + 'px;"></div>';
                if (bHeight > 0)
                    div.innerHTML += tmpl + 'width:' + rPos + 'px; border-right-style:' + (rPos >= 0 ? 'solid' : 'none')
                    + '; border-left-color: darkgreen; border-right-color: darkgreen; opacity: 0.7;'
                    + 'top:' + bTop + 'px; height:' + bHeight + 'px; left:' + lPos + 'px;"></div>';
            }
        });

         var view = this;
         div.onmousemove = function(e) { view.highlightElement([e.pageX , e.pageY]); };
        return div;
    },

    removeGraphTracks: function(idx) {
        var trx = this.m_graphTracks;
        if (!trx) return null;
        if (idx !== undefined) trx = trx.splice(idx, 1);
        else delete this.m_graphTracks;
        Ext.each(trx, function(t) {
            Ext.removeNode(t.tnEl);
            Ext.removeNode(t.div);
            Ext.removeNode(t.lgEl);
            Ext.each(t.labels, function() { Ext.removeNode(this); });
        });
        return trx[0];
    },

    addTracksImage: function(image){
        if (!this.m_Loading) return;
        var graphDiv = Ext.get(this.m_DivId);
        var screen_width = graphDiv.getWidth();
        var szTrackFont = 12;
        var nodes = graphDiv.dom.childNodes;
        var tracks = [];
        var trackConfig = this.m_App.m_Config.TrackConfig;
        var height = 0;
        var siblingName = siblingImg = null;
        for (var i = image.idx - 1; i >= 0; i--) {
            if (this.graphTracksNew[i]) {
                height += this.graphTracksNew[i][0].image.h;
                siblingName = siblingName || this.graphTracksNew[i][0].image.lastNameEl;
                siblingImg = siblingImg || this.graphTracksNew[i][0].image.lastImgEl;
            }
        }
        for (var i in image.areas) { if (!NCBIGBUtils.isNumeric(i)) continue;
            var area = image.areas[i];
            if (!(area.type & (SeqView.AreaFlags.Track | SeqView.AreaFlags.Sequence))) continue;

            var track = {l: area.bounds.l, imgT: area.bounds.t, imgH: image.h - area.bounds.t, divT: area.bounds.t + height, display_name: ''};
            if (area.type & SeqView.AreaFlags.Sequence) {
                track.primers = area.bounds.b - area.bounds.t;
                //skip this area since we have one more
                image.areas.splice(i, 1);
                area = image.areas[i];
            }
            var signature = track.signature = area.signature.split(';')[0];
            track.idx = 0;
            track.image = image;
            Ext.each(trackConfig, function(trk, idx) {
                if (signature != trk.name) return true;
                track.idx = idx;
                track.id = trk.id;
                if (trk.legend) track.legends = trk.legend;
                if (trk.key != 'sequence_track' || trk.check_boxes[0].value) track.display_name = trk.display_name;
                return false;
            });
            trackConfig[track.idx].no_navi = (area.type & SeqView.AreaFlags.NoNavigation);
            this.tmp_showGeneModeButton |= trackConfig[track.idx].category.name == 'Genes';

            if (area.status && area.status.code != 0) track.status = area.status;
            track.area = area;
            area.isSlim = (area.bounds.t == area.bounds.b && track.display_name.length > 0); // slim track (w/o title bar)

            for (var j in this.m_graphTracks) {
                if (this.m_graphTracks[j].id == track.id) {
                    this.removeGraphTracks(j);
                    break;
                }
            }
            if (track.display_name) { // track name div creation
                var newEl = track.tnEl = new Ext.Element(document.createElement('div'));
                var displayNameMod = "";
                Ext.each(trackConfig[track.idx].choice_list, function(i_ch_list) {
                    if (i_ch_list.name == "scale") {
                        if (i_ch_list.curr_value != "linear") {
                            Ext.each(i_ch_list.values, function(i_ch_value) {
                                if (i_ch_value.name == i_ch_list.curr_value) {
                                    displayNameMod = i_ch_value.display_name;
                                    return false;
                                }
                            });
                        }
                        return false;
                    }
                });
                newEl.name = track.display_name;
                if (displayNameMod && displayNameMod.length > 0) {
                    newEl.name += " - " + displayNameMod.toLowerCase() + " scaled";
                }
                if (track.status) {
                    newEl.name += '<div style="color: red; position: relative; top: -5px; font-size:'
                               + (szTrackFont - 2) + 'px; font-family: Arial, san-serif;">Error: '
                               + track.status.msg + '</div>';
                }
                newEl.dom.idx = track.idx;
                newEl.addCls(['sv-drag']);


                if (!siblingName) siblingName = graphDiv.appendChild(newEl);
                else siblingName = newEl.insertAfter(siblingName);

                newEl.applyStyles('text-align: left; position: absolute; top:' + track.divT
                                + 'px; left: 0px; font-family: "Monospace", Courier New, serif; font-size: ' + szTrackFont + 'px;'
                                + (track.area.isSlim ? 'opacity: 0.3; cursor: pointer;' :  ('width:' + this.m_Width +'px;'))
                                + 'overflow: hidden; white-space: nowrap; text-overflow: ellipsis;');
                newEl.insertHtml('afterBegin', newEl.name);
                if (track.area.isSlim) {
                    track.area.bounds.b += newEl.getHeight();
                    track.area.tname = newEl;
                    track.area.bounds.r = newEl.getWidth();
                    newEl.dom.area = track.area;
                    newEl.on('mouseenter',function(e) {
                        this.removeFloatingSelection();
                        if (e.target.area) {
                            
                            if (this.m_tnameDiv) try { this.m_tnameDiv.setStyle('opacity', 0.3); } catch(e) {
                            }
                            this.m_tnameDiv = e.target.area.tname;
                            this.m_Selection = new SeqView.Selection(this, [e.target.area]);
                            this.m_tnameDiv.setStyle('opacity', 1);
                        }
                    }, this);
                } else {
                    newEl.on('mouseover',function(e) {
                        this.highlightElement(e.getXY());
                    }, this);
                }

            }// end of track name div
            if (tracks.length) tracks[tracks.length - 1].imgH -= track.imgH;
            tracks.push(track);

        } //for: image.areas
        for (var i = image.areas.length; i--;) {
            var area = image.areas[i];
            if (!(area.type & SeqView.AreaFlags.Legend)
                || !(area.type & SeqView.AreaFlags.NoSelection)) continue;
            image.areas.splice(i, 1);
            Ext.each(tracks, function(t) {
                if (t.id != area.parent_id) return;
                t.legends = t.legends || [];
                t.legends.push(area);
                return false;
            });
        }
        for (i = 0; i < tracks.length; i++) {
            var track = tracks[i];
            var tname = 'track_' + this.m_Idx + '_' + track.idx;
            if (!track.display_name.length && track.signature == 'sequence') tname += 'SqTrk';
            var wrapDiv = this.createTrackDiv(track, screen_width, tname, siblingImg);
/*            if (!siblingImg) siblingImg = graphDiv.insertFirst(wrapDiv);
            else siblingImg = wrapDiv.insertAfter(siblingImg);*/
            
            track.div = wrapDiv;
            if (track.legends) track.lgEl = wrapDiv.appendChild(this.createLegendsDiv(track, trackConfig));
        }
        this.graphTracksNew[image.idx] = tracks;
        for (var i = image.idx + 1; i < this.graphTracksNew.length; i++) {
            Ext.each(this.graphTracksNew[i], function() {
                this.div.setTop(this.divT += image.h);
                if (this.tnEl) this.tnEl.setTop(this.divT);
            }); 
        }
        image.lastNameEl = siblingName;
        image.lastImgEl = siblingImg;
        this.m_Height += image.data_height;
        if (graphDiv.getHeight() < this.m_Height) {
            graphDiv.setHeight(this.m_Height);
            this.m_View.updateLayout();
        }

        if (this.m_graphTracks && this.m_graphTracks.length && !this.deferGTClean)
            this.deferGTClean = Ext.defer(function () {
                if (this.graphTracksNew) this.removeGraphTracks();
                delete this.deferGTClean;
            }, 1000, this)
    },


    processTracksImages: function(){
        var graphDiv = Ext.get(this.m_DivId);
        document.getElementById(this.m_App.m_DivId).style.cursor = '';
        if (this.deferGTClean) clearTimeout(this.deferGTClean);
        delete this.deferGTClean;
        

        var nodes = graphDiv.dom.childNodes;
        var tracks = [],
            tracks_errors = false;
        var trackConfig = this.m_App.m_Config.TrackConfig;
        var data = this.m_FromCgi;
        var lPos = this.m_Width,
            rPos = 0;

        this.m_Height = 0;

        this.hideGeneModeButton(!this.tmp_showGeneModeButton);
        delete this.tmp_showGeneModeButton;

        this.removeGraphTracks();


        this.m_bottomRulerHeight = this.m_topRulerHeight = data.images[0].areas.shift().bounds.b + 1;

        for (var j = 0; j < data.images.length; j++) {
            var image = data.images[j];
            tracks = tracks.concat(this.graphTracksNew[j]);
            if (image.hairlines && data.images.length > 1) {
                this.m_featMarkers = this.m_featMarkers || {hls: []};
                var arrHLs = image.hairlines;
                arrHLs.sort(function(a, b){return a.pos - b.pos});
                this.m_featMarkers.hls[j] = {hairlines: arrHLs, top: this.m_Height || this.m_topRulerHeight,
                    height: image.h - this.m_topRulerHeight * !j};
                lPos = Math.min(lPos, arrHLs[0].pos);
                rPos = Math.max(rPos, arrHLs[arrHLs.length - 1].pos);
            }
            for (var i = 0; i < image.areas.length; i++) {
                var area = image.areas[i];
                area.bounds.t += this.m_Height;
                area.bounds.b += this.m_Height;
            }
            this.m_Height += image.h;
        }
        this.m_graphTracks = tracks;
        for (var i in tracks) { if (!NCBIGBUtils.isNumeric(i)) continue;
            var track = tracks[i];
            if (track.status) tracks_errors = true;
            if (!track.primers) continue;
            var seq_t = track.area.bounds.t - 5;
            var seq_b = track.area.bounds.t + track.primers - 8;

           //left2right sequence 5' label
           var newEl1 = new Ext.Element(document.createElement('div'));
           var img = '<img data-qtitle="5 \'" data-qwidth="80" data-qtip="5 Prime end" src="' + SeqView.base_url + 'images/5prime_';
           newEl1.insertHtml('afterBegin', img + 'r.gif">');
           newEl1.applyStyles('position:absolute;z-index:10;left:0px;opacity:0.6;top:'
               + ((!this.m_Flip) ? seq_t : seq_b) +'px;');
           //right2left sequence 5' label
           var newEl2 = new Ext.Element(document.createElement('div'));
           newEl2.insertHtml('afterBegin', img + 'l.gif">');
           newEl2.applyStyles('position:absolute;z-index:10;left:' + (graphDiv.getWidth() - 12) + 'px;opacity:0.6;top:'
               + ((!this.m_Flip) ? seq_b : seq_t) +'px;');
           track.labels = [graphDiv.appendChild(newEl1.dom), graphDiv.appendChild(newEl2.dom)];
        } 

        delete this.graphTracksNew;
        if (Ext.isIE8) this.m_graphTracks = [];
        if (this.m_featMarkers) {
            this.m_featMarkers.left = lPos;
            this.m_featMarkers.width = rPos - lPos;
        }

        // adding comment labels for the picture
        for (i = data.areas.length - 1; i >= 0; i--) {
            if (!(data.areas[i].type & SeqView.AreaFlags.Comment)) continue;
            var area = data.areas.splice(i, 1)[0];
            var newEl = new Ext.Element(document.createElement('div'));
            newEl.insertHtml('afterBegin', area.label);
            var styles = 'opacity:0.6;position:absolute;font-family: "Monospace", Courier New, serif; font-size: 12px;';//should always be left:0px;
            if (area.bounds.t) styles += 'top:' + area.bounds.t + 'px;';
            if (area.bounds.l == -1) styles += 'left:' + area.bounds.r + 'px;';
            if (area.bounds.l == 1) styles += 'right:' + (el.getRight() - area.bounds.r) + 'px;';
            if (area.type & SeqView.AreaFlags.DrawBackground) {
                styles += 'background-color:rgb('
                    + (this.m_App.m_Config.Options.curr_color == 'Color' ? '180,180,240' : '195,195,195') + ');';
            }
            Ext.each(tracks, function(t) {
                if (t.id != area.parent_id) return;
                t.labels = t.labels || [];
                t.labels.push(graphDiv.appendChild(newEl.dom));
                newEl.applyStyles(styles);
                return false;
            });           

            if (area.tooltip) var tip = new Ext.ToolTip({target: newEl, html: area.tooltip});
            newEl.on('mouseover', function(e) { this.highlightElement(e.getXY()); }, this);
        }
        if (tracks_errors) this.m_App.fireEvent('tracks_errors', tracks);
        if (!this.m_Locator && this.m_App.m_Panorama) {
            this.m_Locator = new SeqView.Locator(this, this.m_Color, true);
            this.m_Locator.setColor(this.m_Color);
        }
    },

    applyResultsOnLoad: function(options){
        this.m_App.watchdogStop();
        if (!this.m_Loading) return;
        this.m_Loading = false;

        var ts = new Date().getTime();
        var appTS = this.m_App.m_timeStamps;
        var stat = {'SV_graphics_time':  ts - appTS.loadImagesTS}
        if (appTS.prior2seqconfig) {
            stat['sv-event'] = 'initialization';
            stat['SV_config_time'] = appTS.loadImagesTS - appTS.prior2seqconfig;
            stat['SV_config_graphics_time'] = ts - appTS.prior2seqconfig;
            delete appTS.prior2seqconfig;
        }
        this.ping(stat);

        var load_cmp = Ext.getCmp('seq-view-loading-' + this.m_Idx);
        if (load_cmp) load_cmp.enable();

        Ext.get(this.m_DivId).unmask();

        var the_div = Ext.get(this.m_DivId);
        the_div.setStyle('background', '');

        var screen_width = the_div.getWidth(); //this.getScreenWidth();

        var chNodes = the_div.dom.childNodes;
        var trackConfig = this.m_App.m_Config.TrackConfig;
        var gview = this;
        for (var i = chNodes.length - 1; i >= 0; i--) {
            var node = chNodes[i];
            if (typeof node.idx == 'undefined' || !node.textContent || !node.textContent.length) continue;
            if (node.style.opacity || this.m_App.m_PermConfId || Ext.isIE8) continue;
            node.style.width = screen_width + 'px';
            var tpl;
            var no_navi = trackConfig[node.idx].no_navi;
            if (!no_navi) { // to avoid overlapping jump buttons
                tpl = new Ext.Template ('<div style="width:' + ((screen_width >> 1) - 26)
                    + 'px;overflow: hidden; white-space: nowrap; text-overflow: ellipsis;">' + node.innerHTML + '</div>');
                node.textContent = '';
                tpl.append(node, {}, true).dom;
            }
            if (!this.m_App.m_NoConfDlg) { // X-button adding
                tpl = new Ext.Template('<button style="height:13px; width:13px; position: absolute; top:-1px; left:'
                    + (screen_width - 13) + 'px; cursor: pointer; background-image: url('+ SeqView.base_url
                    + 'images/hide_track_red.gif); padding: 0px;'
                    + 'border-radius: 3px; border: 0px; background-color: transparent;'
                    + 'background-position:center;background-repeat:no-repeat;z-index:20;" title="Hide track"></button>');
                var butt = tpl.append(node, {}, true).dom;
                butt.onmouseover = function(e) { this.style.border = ''; };
                butt.onmouseout = function(e) { this.style.border = '0px'; };
                butt.ontouchstart = butt.onclick = function(e) {
                    e.stopPropagation();
                    e.preventDefault();
                    gview.pingClick('2-0-1');
                    trackConfig[this.parentNode.idx].shown = false;
                    if (this.m_SelectedSig && !this.m_featMarkers)
                        SeqView.TM.Common.updateSeqViewApp(SeqView.TM.processTracksInfo(trackConfig), gview.m_App);
                    else 
                        gview.hideTrack(this.parentNode.idx);
                };
            }
            if (no_navi) continue;

            tpl = new Ext.Template('<button style="height:11px; position: absolute; top:2px; left:'
                + ((screen_width >> 1) - 26) + 'px; cursor: pointer; padding: 0px; border: 0px;'
                + 'border-radius: 1px; z-index:20;" title="Jump left"></button>');
            var lButt = tpl.append(node, {}, true).dom;
            tpl = new Ext.Template('<button style="height:11px; position: absolute; top:2px; left:'
                + ((screen_width >> 1) + 8) + 'px; cursor: pointer; padding: 0px; border: 0px;'
                + 'border-radius: 1px; z-index:20;" title="Jump right"></button>');
            var rButt = tpl.append(node, {}, true).dom;;
            lButt.posX = 0; rButt.posX = -18;
            lButt.onmouseover = rButt.onmouseover = function(e) {
                if (gview.tmp_gotoFeat || gview.m_Loading) return;
                this.style.opacity = 0.9;
            }
            lButt.onmouseout = rButt.onmouseout = function(e) {
                if (gview.tmp_gotoFeat || gview.m_Loading) return;
                this.style.width = '18px';
                this.style.opacity = 0.5;
                this.style.background = 'url('+ SeqView.base_url + 'images/jumps.png) ' + (this.posX + 'px 0px');
            }
            lButt.onmouseout();
            rButt.onmouseout();
            lButt.up = !(rButt.up = !gview.m_Flip);
            lButt.ontouchstart = rButt.ontouchstart = lButt.onclick = rButt.onclick = function(e) {
                e.stopPropagation();
                e.preventDefault();
                gview.pingClick('2-0-3');
                if (gview.tmp_gotoFeat || gview.m_Loading) return;
                this.style.background = SeqView.ExtIconLoading + ' 0px 0px';
                this.style.width = '11px';
                this.style.backgroundSize = 'cover';
                gview.gotoFeature(this.parentNode.idx, (this.up ? 'next' : 'prev') + (e.ctrlKey || e.metaKey ? '_group' : ''),
                    function(){ e.target.onmouseout(); });
            }
        }
        var trx = this.m_graphTracks;

        var id = 'tRuler' + this.m_Idx;
        var wrapDiv = Ext.get(id) || new Ext.Element(document.createElement('div'));
        wrapDiv.set({id: id, style: 'position:absolute;overflow:hidden;top:0;width:'
            + this.m_Width + 'px;height:' + this.m_topRulerHeight + 'px;left:' + this.m_ScrollPix + 'px;'});
        var img = Ext.get(id + 'img') || new Ext.Element(document.createElement('img'));
        img.set({id: id + 'img', draggable: false, src: this.m_FromCgi.images[0].img_url});
        the_div.appendChild(wrapDiv);
        wrapDiv.appendChild(img);
        if (this.m_bottomRulerHeight >= 0) {
            id += 'B';
            var wrapDiv = Ext.get(id) || new Ext.Element(document.createElement('div'));
            wrapDiv.set({id: id, style: 'position:absolute;top:' + this.m_Height + 'px; width:'
                + this.m_Width + 'px; height:' + this.m_bottomRulerHeight + 'px;left:' + this.m_ScrollPix + 'px;'});
            var img = Ext.get(id + 'img') || new Ext.Element(document.createElement('img'));
            img.set({id: id + 'img', draggable: false, src: this.m_FromCgi.images[0].img_url});
            this.m_Height += this.m_bottomRulerHeight;
            the_div.appendChild(wrapDiv);
            wrapDiv.appendChild(img.dom);
        }
        the_div.setHeight(this.m_Height);
        this.m_vLine.setHeight(this.m_Height);

        if (this.m_FromCgi.areas && this.m_FromCgi.areas.length > 0) {
            this.m_App.addCustomFeatureFlags({view: this, areas: this.m_FromCgi.areas, callback: this.highlightFlaggedFeatures});
        }
        if (this.m_featMarkers)
            the_div.appendChild(this.drawFeatMarkers());
        if (!options || options.update_locator !== false) this.m_App.updateLocator(this);
        this.m_View.updateLayout();
        this.pushToHistory();

        this.m_App.updateMarkersPos(this);
        this.m_App.updateReflections();

        if (this.m_SelectedRangeSet) {
            this.selectRangeSet(this.m_SelectedRangeSet);
            this.m_SelectedRangeSet = null;
        }

        this.blinkAreas(options);

        if (!this.m_App.m_DialogShown) this.m_App.resizeIFrame();
        this.m_App.notifyViewLoaded(this);
        this.m_App.updateMarkersSize(this);
    },

    createLegendsDiv: function(track, trackConfig) {
        var legEl = new Ext.Element(document.createElement('div'));
        var bLen = 30;
        var qview_width = this.getScreenWidth();
        var leftOffset = (this.m_Width - qview_width) >> 1;
        var view = this;
        var img = track.image,
            minT = img.h,
            maxB = 0;
        Ext.each(track.legends, function(lg) {
            if (lg.bounds) {
                lg.label = lg.signature;
                lg.color = 'rgba(' + lg.color.replace(/ /g, ',') + ')';
                lg.opacity = 100;
                minT = Math.min(minT, lg.bounds.t);
                maxB = Math.max(maxB, lg.bounds.b);
            }
            else {
                lg.opacity = parseInt(lg.opacity);
                Ext.each(img.areas, function() {
                    if (!((this.type & SeqView.AreaFlags.Legend)
                        && lg.id == this.id && this.parent_id == trackConfig[track.idx].id)) return true;
                    lg.bounds = this.bounds;
                    minT = Math.min(minT, this.bounds.t);
                    maxB = Math.max(maxB, this.bounds.b);
                });
            }
        });
        legEl.applyStyles('position: absolute; left: 0px;'
                        + 'width:' + qview_width +'px; height:' + (maxB - minT) + 'px;'
                        + 'top:' + (track.imgH - maxB + minT - 3) + 'px;');
        Ext.each(track.legends, function(lg, idx) {
            if (!lg.bounds) return true;
            var ib = lg.bounds;
            var tmpl = new Ext.Template('<div class="sv-drag" style="left:' + (ib.l - leftOffset) + 'px; top:' + (ib.t - minT)+ 'px;'
                + 'width:' + (ib.r - ib.l) + 'px; height:' + (ib.b - ib.t) +'px; position: absolute;'
                + (lg.parent_id ? '' : 'cursor: pointer') + ';">'
                + '<div style="width:' + bLen + 'px; height: 8px; position:relative;'
                + 'background-color:' + lg.color + '; opacity:' + (lg.opacity/100)+ ';"></div>'
                + '<div style="color:blue; font-size: 10px; top: -10px; left:' + (bLen + 5)
                + 'px; width:' + (ib.r - ib.l - bLen - 5) + 'px; overflow: hidden; white-space: nowrap; text-overflow: ellipsis; position:relative;">'
                + lg.label + '</div></div>');
            var legTrk = tmpl.append(legEl, {}, true).dom;
            if (lg.parent_id) return;
            legTrk.onclick = function(e) {
                e.stopPropagation();
                e.preventDefault();
                view.pingClick('4-1', 'Graph_Overlay');
                SeqView.TM.modifyLegendSettings(view, track, idx);
            }
        });
        return legEl.dom;
    },

    createTrackDiv: function(track, width, id, sibling) {
        var wrapDiv = new Ext.Element(document.createElement('div'));
        wrapDiv.set({id: id});
        wrapDiv.applyStyles('top: ' + track.divT + 'px; width: ' + width
                          + 'px; height: ' + track.imgH + 'px; position:absolute; overflow: hidden;');
        if (sibling) wrapDiv.insertAfter(sibling);
        else sibling = Ext.get(this.m_DivId).insertFirst(wrapDiv);
/*        
        var cnvs = wrapDiv.insertFirst(document.createElement('canvas'));
        cnvs.height = track.imgH;
        cnvs.width = track.img
        var cntx = cnvs.dom.getContext('2d');
        */
        var img = document.createElement('img');
        img.id = 'img' + id;
        img.setAttribute('class', 'sv-drag sv-highlight sv-dblclick');
        img.setAttribute('style', 'position:relative; overflow:hidden; visibility:visible;'// + 'display:none;';
           + 'left:' + this.m_ScrollPix + 'px;top:' + (-track.imgT) + 'px;');
        img.src = track.image.img_url;
        wrapDiv.insertFirst(img);
//        var cntx = wrapDiv.dom.getContext('2d');
//        cntx.drawImage(img,  this.m_ScrollPix, -track.imgT);
        return wrapDiv;
    },

    loadFailure: function(text, res){
        if (res) this.m_App.watchdogReport(res);
        this.m_View.setTitle('Image loading error: ' + text);
        this.m_Loading = false;
        this.mf_HistoricalBack = false;
    },

    gotoFeature: function(track_idx, dir, callback) {
        callback = callback || function(){};
        this.tmp_gotoFeat = true;
        var vLen = this.m_VisLenSeq,
            vFrom = this.m_VisFromSeq,
            lenSeq = this.m_LenSeq,
            fromSeq = this.m_FromSeq;
        var curMarker = this.m_App.getMarkersInfo().findMarkerByName(SeqView.MarkerNav);
        var sPos = (curMarker) ? curMarker.seq_pos : (vFrom + (vLen >> 1));

        var processError = function(data, text, res) {
            delete this.tmp_gotoFeat;
            console.log('gotoFeature error:' + text);
            callback();
       }

        var processResponse = function(data, text, res) {
            if (data.job_status) {
                if (data.job_status == 'failed' || data.job_status == 'canceled') {
                    processError(null, data.job_status);
                } else {
                    Ext.defer(SeqView.App.simpleAjaxRequest, 500, this, [{
                        url: this.m_App.m_CGIs.Graphic + '?job_key=' + data.job_id,
                        context: this,
                        success: processResponse,
                        error: processError}]);
                }
                return;
            }
            if (typeof data.new_pos !== 'undefined') {
                var dPos = data.new_pos;
                var bOpt = {blinkData: {track_idx: track_idx, data: data}};
                var newVisFrom = dPos - (vLen >> 1);
                if (newVisFrom < 0) newVisFrom = 0;
                if ((vLen < this.m_App.m_SeqLength) && (dPos < vFrom || dPos > (vFrom + vLen))) {
                    if (newVisFrom < fromSeq || (fromSeq + lenSeq) < (newVisFrom + vLen)) {
                        this.gotoPosRange([dPos, 0], true, bOpt);
                    } else {
                        var shift = Math.round((vFrom - newVisFrom) * this.m_bpPix);
                        if (!this.scrollViewTo(this.m_ScrollPix + shift, shift < 0 ? SeqView.PAN_RIGHT : SeqView.PAN_LEFT, false, bOpt))
                            this.blinkAreas(bOpt);
                    }
                } else {
                    this.blinkAreas(bOpt);
                }
                if (curMarker) {
                    curMarker.setSeqPos(dPos, false);
                    curMarker.m_ToolTip.updateMarkerToolTipContent(null, curMarker.span, dPos);
                }
                else this.m_App.getMarkersInfo().addMarker([dPos], SeqView.MarkerNav, false, 'navy');
            }
            delete this.tmp_gotoFeat;
            callback();
        }
        var data = {track: SeqView.TM.tracksArrayToString(this.m_App.m_Config.TrackConfig[track_idx], true, true),
                   pos: sPos, select: this.m_SelectedSig, dir: dir, id: this.m_App.GI, assm_context: this.m_App.m_AssmContext,
                   navi: 1, from: fromSeq, len: lenSeq, width: this.m_Width};
        Ext.apply(data, this.m_App.m_Config.visualOptionsUrl());
        this.m_App.AjaxRequest({ url: this.m_App.m_CGIs.Graphic, data: data, context: this,
            success: processResponse,
            error: processError});
    },

    blinkAreas: function(options) {
        if (!options || !options.blinkData) return;
        var data = options.blinkData.data;
        if (!data.bounds || !data.bounds.length) return;
        var track, idx = options.blinkData.track_idx;
        Ext.each(this.m_graphTracks, function(t) { track = t; return t.idx != idx; });
        var offset = this.m_FromSeq + (this.m_Flip ? this.m_LenSeq : 0);
        track.blinkData = data;
        Ext.each(data.bounds, function(db) {
            db.l = Math.abs(Math.round((db.l - offset) * this.m_bpPix));
            db.r = Math.abs(Math.round((db.r - offset) * this.m_bpPix));
            db.t += track.divT + 1;
            db.b += track.divT + 1;
            tpl = new Ext.Template('<div style="position:absolute; border: 1px solid red;"></div>');
            db.elem = SeqView.createSelectionElement(this, {bounds: db}, tpl);
        }, this);
        function blink(bnum) {
            var stop = bnum-- < 0 || track.blinkData != data;
            if (!stop) Ext.defer(blink, 300, this, [bnum]);
            Ext.each(data.bounds, function(r) {
                if (!stop) r.elem.setVisible(bnum % 2);
                else r.elem.destroy();
            });
        }
        blink(3);
     },

    highlightFlaggedFeatures: function() {
        if (this.m_FromCgi.areas) {
            this.removeFeatMarks();
            Ext.each(this.m_FromCgi.areas, function(a) {
                if( (a.type & SeqView.AreaFlags.Dirty) !== 0) {
                   var sh = new SeqView.SelectionHighlighter(this, a, true);
                   if (!this.m_HighLightedFeature) this.m_HighLightedFeature = [];
                   this.m_HighLightedFeature.push(sh);
                }
            },this);
        }
    },

    shiftGraphTrack: function(trk, shift) {
        shift = shift || trk.shift || 0;
        trk.div.setTop(trk.divT += shift);
        trk.area.bounds.t += shift;
        trk.area.bounds.b += shift;
        if (trk.tnEl) trk.tnEl.setTop(trk.divT);
        Ext.each(trk.labels, function() { this.setTop(this.getTop(true) + trk.shift); });
    },

    hideTrack: function(cfgIdx) {
        var gtIdx;
        this.m_graphTracks.some(function(t, i) { gtIdx = i; return t.idx == cfgIdx; });
        this.m_App.m_actTracks = '';
        var trk = this.removeGraphTracks(gtIdx);
        if (trk.image.hairlines && this.m_featMarkers){
            this.m_featMarkers.div.parentElement.removeChild(this.m_featMarkers.div);
            delete this.m_featMarkers;
        };
        this.m_App.fireEvent('configuration_changed', this.m_App);
        if (this.m_featMarkers) {
            this.refresh();
            return;
        }
        
        Ext.each(this.m_graphTracks, function(t, i) { if (i >= gtIdx) this.shiftGraphTrack(t, -trk.imgH); }, this);

        Ext.each(this.m_FromCgi.areas, function(area, idx) {
            if (area == trk.area || area.parent_id == trk.id) { 
                this.splice(idx, 1);
                return;
            }
            if (!area.parent_id || area.bounds.t < trk.divT) return;
            area.bounds.t -= trk.imgH;
            area.bounds.b -= trk.imgH;
        }, this.m_FromCgi.areas, true);

        var graphDiv = Ext.get(this.m_DivId);
        graphDiv.setHeight(this.m_Height -= trk.imgH);
        var ruler = Ext.get('tRuler' + this.m_Idx + 'B');
        if (ruler) ruler.setTop(this.m_Height - this.m_bottomRulerHeight);
        this.selectRangeSet(this.getSelectedRangeSet());
        this.m_App.updateMarkersSize(this);
        this.updateStatusBar('tracks', 'Tracks shown: ' + + this.m_graphTracks.length + '/' + this.m_App.m_Config.TrackConfig.length)
        this.m_View.updateLayout();
    },


    shiftDownTrack: function(trkIdx) {
        this.m_App.m_actTracks = '';
        var trx = this.m_graphTracks;
        var src = trx[trkIdx], dst = trx[trkIdx + 1];
        dst.shift = -src.imgH;
        src.shift = dst.imgH;

        this.shiftGraphTrack(src);
        this.shiftGraphTrack(dst);
        
        trx.splice(trkIdx + 2 , 0, trx[trkIdx]);
        trx.splice(trkIdx, 1);

        Ext.each(this.m_FromCgi.areas, function(area) {
            var shift = dst.shift;
            if (area.parent_id == src.id) shift = src.shift;
            else 
                if (area.parent_id != dst.id) return true;
            area.bounds.t += shift;
            area.bounds.b += shift;
        }, this);
    },

    clearCache: function() {
        this.m_PrevCgi= null; // clear next/prev cache
        this.m_NextCgi = null;
        this.m_UrlFrom = null;
        this.m_UrlTo = null;
    },

    scrollElements: function(delta) {
        this.m_App.scrollReflections(this, delta);
        this.m_App.scrollMarkers(this, delta);
        if (this.m_featMarkers) {
            this.m_featMarkers.div.style.left = this.m_featMarkers.left + this.m_ScrollPix + 'px';
        }

        if (this.m_HighLightedFeature) {
            Ext.each(this.m_HighLightedFeature,function(sf) {
                sf.movePix(delta);
            });
        }
        if (this.m_Selection) { this.m_Selection.movePix(delta); }
        var sel_set = this.m_RangeSelectionSet;
        if (sel_set && sel_set.length) {
            var l = sel_set.length;
            for (var x = 0; x < l; ++x)
                sel_set[x].movePix(delta);
        }
    },

    updateImagePosition: function(new_pos){
        var delta = this.m_ScrollPix - new_pos;
        if (delta == 0) return false;

        delete this.m_ReqLen;

        this.m_ScrollPix = new_pos;
        this.scrollElements(delta);

        if (!this.m_Loading) {//ImageShown) {
            var the_div = Ext.get(this.m_DivId);
            var npos = this.m_ScrollPix + the_div.getX();
            Ext.each(the_div.dom.children, function () {
                 if (this.id.indexOf('track_') == 0)
                     Ext.get(this.children[0]).setX(npos);
            });
        }
        return true;
    },

//////////////////////////////////////////////////////////////////////////
// scrollView:

    scrollViewTo: function(new_pos, dir, drag, options){
// returns TRUE if a new image was requested to load
        drag = drag || false;

        this.m_UrlFrom = null;
        this.m_UrlTo = null;

        var screen_width = this.getScreenWidth();

        var new_img_pos = Math.max(Math.min(0, new_pos), screen_width-this.m_Width);

        if (!this.updateImagePosition(new_img_pos)) return false;

        var vis_range = this.toSeq();
        var preloadMargin = Math.min(SeqView.ChunkWidth, this.m_Width/10);

        if (!this.m_Loading &&
             (  (new_pos < screen_width - this.m_Width + preloadMargin && dir == SeqView.PAN_RIGHT)
             || (new_pos > -preloadMargin                              && dir == SeqView.PAN_LEFT)))
        {
            var new_from, new_len;
            var bpp = Math.max(SeqView.MinBpp, this.m_bpPix);

            var loading_needed = true;
            var moving_to_seq_end = this.getFlip() ? (dir == SeqView.PAN_LEFT) : (dir == SeqView.PAN_RIGHT);
            if (moving_to_seq_end) {
                if (this.m_FromSeq + this.m_LenSeq == this.m_App.m_SeqLength - 1) {
                    loading_needed = false;
                }            
            } else {
                if (this.m_FromSeq == 0) { loading_needed = false; } // start of sequence reached - no loading
            }
            if (loading_needed) {
                new_from = vis_range[0];// + vis_range[2] / 2;
                new_len = vis_range[2];
                this.m_Loading = true;
                this.m_ScrollPixBeforeLoad = this.m_ScrollPix; // save scroll offset at the time of load request
                var opt = options || {};
                opt.from_ui = opt.transitional = true;
                this.loadImage(new_from, new_len, opt);
                return true;
            }
        }
        if (!drag) {
            this.updateTitle({from_ui: true, transitional: drag }); // and save visible range
            this.pushToHistory(vis_range);
        }
        return false;
    },


    updateTitle: function(options) {
        var range = this.toSeq();
        var from = (range[0] + 1);
        var to = (range[1] + 1);

        if (this.getFlip()) {
            var t = from;
            from = to;
            to = t;
        }

        var units = (this.m_App.m_ViewParams['acc_type'] == 'protein' ? 'residues' : 'bases')
        var title = this.m_App.m_ViewParams['id'] + ':&nbsp;';
        title += from.shorten() + '..' + to.shorten() + '&nbsp;(' + range[2].shorten();
        title += (units.charAt(0) == 'b' ? 'bp' : 'r') + ')';
        if (this.canFlip() && this.getFlip()) title += ' C';

        var tbar = this.m_View.getDockedItems('toolbar[dock="top"]')[0];
        if (tbar) {
            var tbtitle = tbar.down('#tbtitle');
            if (tbtitle) {
                var tiptitle = '<b>' + this.m_App.m_ViewParams['id'];
                tiptitle += ': ' + this.m_App.m_Config.SeqInfo.title + '</b><br><br>';
                tiptitle += from.commify() + '&nbsp;-&nbsp;' + to.commify() + '<br>';
                tiptitle += range[2].commify() + '&nbsp;' + units + '&nbsp;shown';
                if (this.canFlip())
                   tiptitle += ',&nbsp;' + (this.getFlip() ? 'reverse complement' : 'positive strand');

                tbtitle.setTooltip(tiptitle);

                if (this.getScreenWidth() < 650){
                    tbtitle.setWidth(undefined);
                    tbtitle.setText('');
                    tbtitle.setIconCls('xsv-seq-logo');
                } else {
                    tbtitle.setText('<b style="font-size:x-small;">' + title + '</b>');
                    tbtitle.setIconCls('');
                }
            }

            //taking care of showing either "Find on Sequence" button with text field or just button
            var searchButton = tbar.down('#btnSearch');
            if (searchButton) searchButton.hide();
            var tbgts = tbar.down('#gotoBox');
            var tbgtsb = tbar.down('#btnFind');
            if (tbgts && tbgtsb) {
                tbgts.show();
                tbgtsb.show();
            } else {
                if (searchButton) searchButton.show();
            }
            var tbtoolsb = tbar.down('#btnTools');
            if (tbtoolsb) tbtoolsb.setText('Tools');
            var tools = tbar.items.items;
            if ((tbgts && tbgtsb) && tbar.el.getRight() < tools[tools.length - 1].el.getRight()) {
                tbgts.hide();
                tbgtsb.hide();
                if (searchButton) searchButton.show();
            }
            if (tbtoolsb && tbar.el.getRight() < tools[tools.length - 1].el.getRight()) {
                tbtoolsb.setText(''); // Tools menu
            }
        }

        this.m_VisFromSeq = range[0];
        this.m_VisLenSeq = range[2];
/*
        if (options && options.from_ui){
            this.m_App.fireEvent('ui_visible_range_changed', this);
        } else {
            this.m_App.fireEvent('api_visible_range_changed', this);
        }
*/
        var flags = {};
        flags.from_ui = options && options.from_ui;
        flags.transitional = options && options.transitional;
        flags.fire = options && options.fire;

        this.m_App.fireEvent( 'visible_range_changed', this, flags );
    },


    pushToHistory: function(range) {
        if( this.m_InDragAction ) return;
        range = range || this.toSeq(); 
        var flip = this.getFlip();
        var MAX_HISTORY_LENGTH = 30;
        var last = this.m_History[this.m_History.length - 1];

        if (this.mf_HistoricalBack || (last && Math.abs(range[0] - last.from) <= 1 && range[2] == last.len && flip == last.flip) ){
            this.mf_HistoricalBack = false;
            return;
        }
        if (this.m_History.push({from: range[0], len: range[2], flip: flip}) > MAX_HISTORY_LENGTH)
            this.m_History.shift();

        this.mf_HistoricalBack = false;
        this.updateHistoryButtons();
    },

    stepHistory: function(options ){
        this.m_History.pop();
        var prev = this.m_History[this.m_History.length - 1];
        if (prev) {
            this.mf_HistoricalBack = true;
            this.setFlip(prev.flip);
            this.startImageLoading(prev.from, prev.len, options);
        }
        this.updateHistoryButtons();
    },

    updateHistoryButtons: function() {
        var btn_prev = this.m_View.down('#histPrev');
        if (btn_prev) {
            var prev = this.m_History[this.m_History.length - 2];
            btn_prev[prev ? 'enable' : 'disable']();
            var title = 'Back';
            if (prev) {
                var to = prev.from + 1 + (prev.flip ? 0: prev.len);
                var from = prev.from + 1 + (prev.flip ? prev.len : 0);  
                title += ' to ' + from.commify() + '&nbsp;-&nbsp;' + to.commify()
                    + '&nbsp;(' + prev.len.commify() + '&nbsp;'
                    + (this.m_App.m_ViewParams['acc_type']=='protein' ? 'residues' : 'bases');
                if (this.canFlip()) title += ',&nbsp;' + (prev.flip ? 'negative' : 'positive');
                title += ")";
            }
            btn_prev.setTooltip(title);
        }
    },

    moveTo: function(vis_from, vis_len, opts){
        if (vis_len == this.m_VisLenSeq) {
            var cur_vis_from_pix = Math.round(this.seq2Pix(this.m_VisFromSeq));
            var new_vis_from_pix = Math.round(this.seq2Pix(vis_from));
            var delta_x = cur_vis_from_pix - new_vis_from_pix;
            if (delta_x == 0) return;

            var cur_x = Math.round(this.m_ScrollPix);
            cur_x = cur_x + delta_x;
            var rest_len = this.m_Width + this.m_ScrollPix - this.getScreenWidth();
            if ( cur_x > 0 || (delta_x < 0 && -delta_x > rest_len) ){
                this.startImageLoading(vis_from, vis_len, opts);
            } else {
                // No need to check as it is done above
                this.updateImagePosition(cur_x);
                this.updateTitle(opts);
            }
        } else {
            this.startImageLoading(vis_from, vis_len, opts);
        }
    },

    syncToLocator: function() {

        if( !this.m_Locator ) return;

        var vis_from = this.m_App.m_Panorama.toSeq( this.m_Locator.getLeft( true ) );
        var vis_to;

        if( this.m_Locator.m_Action == 'Drag' ){
            vis_to = vis_from + this.m_VisLenSeq - 1; // keep length

        } else {
            vis_to = this.m_App.m_Panorama.toSeq(this.m_Locator.getLeft(true) + this.m_Locator.getWidth());
        }

        var len = vis_to - vis_from +1;

        this.moveTo( vis_from, len, { update_locator: false, from_ui: true } );
    },


    checkLocatorWidth: function(width) {
        var len = this.m_App.m_Panorama.toSeq(width+3);
        return SeqView.MinBpp < len/this.getScreenWidth();
    },


    hitTest: function(page_xy) {
        if (this.m_Loading || this.m_InDragAction)
            return null;

        var the_area = [];
        var elem_xy = Ext.get(this.m_DivId).getXY();
        var xx = page_xy[0] - elem_xy[0];
        var slim = this.m_slimMode;
        var scrollPix = this.m_ScrollPix;
        var yy = page_xy[1] - elem_xy[1] + 3;// - config['top_offset'];
        var isFlip = this.getFlip();
        Ext.each(this.m_FromCgi.areas, function(area) {
            if (!area.signature) return;
            var bounds = area.bounds;
            var x = xx - ((slim && (area.type & (SeqView.AreaFlags.Track | SeqView.AreaFlags.Sequence)) > 0) ? 0 : scrollPix);
            var left  =  isFlip ? bounds.r : bounds.l;
            var right = isFlip ? bounds.l : bounds.r;
            if (x >= left - 1 && x <= right + 1 && yy >= bounds.t && yy < bounds.b) {
                the_area.push(area);
            }
        });

        return the_area.length ? the_area : null;
    },

    highlightElement: function(page_xy) {
        var areas = this.hitTest(page_xy);
//console.log(page_xy[1], (areas ? areas[0].bounds : ''));
        if (this.m_slimMode && areas && (areas[0].type  & (SeqView.AreaFlags.Track | SeqView.AreaFlags.Sequence))) {
            areas.shift();
            if (!areas.length) return;
            if (this.m_Selection && areas[0].signature != this.m_Selection.areas[0].signature)
                this.removeFloatingSelection();
        }
        if (areas && (!this.m_Selection || this.m_Selection.removed)) {
            this.m_Selection = new SeqView.Selection(this, areas);
        } else {
            if (!this.m_Selection || !areas || areas[0].signature != this.m_Selection.areas[0].signature)
                this.removeFloatingSelection();
        }
    },

    onMouseDown: function(e) {
        if (this.m_ContextMenu) this.m_ContextMenu.destroy();
        if ((e.type == 'mousedown' && e.button) || this.m_Loading) return;
        this.m_XY = e.getXY();

        if (e.target.id.indexOf('Ruler') == 1) {
            // Add range selection
            e.stopPropagation();
            if (!e.ctrlKey) this.removeSelectionsWithNoPinnedToolTips();
            this.m_InDragAction = true;
            new SeqView.RangeSelection(this, [this.m_XY[0]], e);
            this.m_XY = null;
            return;
        }
        if (e.type != 'mousedown')
            this.m_deferredContext = Ext.defer(this.showContextMenu, 2000, this, [this.m_XY]);

        if (Ext.fly(e.getTarget()).hasCls('sv-drag')) {
            this.m_InDragAction = true;
            e.stopEvent();
            Ext.fly(this.m_DivId).setStyle('cursor', 'move');

            var view = this;
            var onMove = function(e) {view.onMouseMove(new Ext.event.Event(e ? e : window.event));}
            var onEnd = function(e) {view.onMouseUp(new Ext.event.Event(e ? e : window.event));}
            if (e.button == 0) {
                this.m_DocMouseMove = document.onmousemove;
                this.m_DocMouseUp = document.onmouseup;
                document.onmousemove = onMove;
                document.onmouseup = onEnd;
            } else {
                this.m_DocTouchMove = document.ontouchmove;
                this.m_DocTouchUp = document.ontouchend;
                document.ontouchmove = onMove;
                document.ontouchend = onEnd;
            }
        }
    },

    onMouseUp: function(e) {
        this.m_InDragAction = false;
        if (!this.m_XY) return;
        if (this.m_deferredContext) {
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        Ext.fly(this.m_DivId).setStyle('cursor', '');//'default');
        e.stopEvent();

        if (this.m_targetTrackClone) {
            Ext.removeNode(this.m_targetTrackClone.dom);
            var app = this.m_App;
            delete this.m_targetTrackClone;
            delete app.m_currentTrackSetId;
            app.fireEvent('configuration_changed', app);
            this.pingClick('2-0-0');
            if (this.m_SelectedSig) {
                this.refresh();
            }
        }
        if (this.m_scrolled) {
            this.updateTitle({from_ui: true, transitional: false, fire: true });
            delete this.m_scrolled;
            clearTimeout(SeqView.scrolled);
            SeqView.scrolled = setTimeout(function() { delete SeqView.scrolled; }, 500);
            if (this.m_ReqFrom != this.m_VisFromSeq) this.pingClick('2-0-2');
        }
        this.m_XY = null;
        if (e.button == 0) {
            document.onmousemove = this.m_DocMouseMove;
            document.onmouseup = this.m_DocMouseUp;
        } else {
            document.ontouchmove = this.m_DocTouchMove;
            document.ontouchend = this.m_DocTouchUp;
        }
        this.m_DocMouseMove = this.m_DocMouseUp = this.m_DocTouchMove = this.m_DocTouchUp = null;

        this.pushToHistory();
    },

    onMouseMove: function(e) {
        NCBIGBUtils.ClearBrowserSelection();

        if (this.m_InDragAction) {
            if (!this.m_XY) { return; }
            if (e.button > 0) {
                this.onMouseUp(e)
                return;
            }

            var sensivity = 2;
            var xy = e.getXY();
            var delta_x = this.m_XY[0] - xy[0];
            var delta_y = this.m_XY[1] - xy[1];
            // ignore vertical moving if there is a permanent ControlPanel
            if (this.m_App.m_PermConfId) delta_y = 0;
            if (delta_x == 0 && delta_y == 0) return;
            var XminusY = Math.abs(delta_x) - Math.abs(delta_y);
            if (Math.abs(XminusY) <= sensivity) return;
            if (this.m_deferredContext) {
                clearTimeout(this.m_deferredContext);
                this.m_deferredContext = 0;
            }

            if (XminusY > 0 ) {
                var cur_x = this.m_ScrollPix;
                cur_x = cur_x - delta_x;
                this.m_XY = xy; // save new values
                this.scrollViewTo(cur_x, delta_x > 0 ? SeqView.PAN_RIGHT : SeqView.PAN_LEFT, true);
                this.m_scrolled = true;
            } else { // tracks moving processing
                var tracks = this.m_graphTracks;
                var the_div = Ext.get(this.m_DivId);
                var src = -1, dst;
                var pY = this.m_XY[1] - the_div.getY();

                for (var i = 0; i < tracks.length; i++) {
                    if (pY > tracks[i].divT && pY < tracks[i].divT + tracks[i].imgH) {
                        src = i;
                        break;
                    }
                }

                if (src >= 0) {
                    var srcY = tracks[src].divT - delta_y;
                    if (delta_y > 0) {// moving up
                        for (dst = 0; dst < src; dst++) {
                            pY = tracks[dst].divT - srcY;
                            if (Math.abs(pY) <= sensivity) break;
                        }
                    } else {// moving down
                        srcY = tracks[src].divT+ tracks[src].imgH - delta_y;
                        for (dst = tracks.length - 1; dst > src; dst--) {
                            pY = tracks[dst].divT+ tracks[dst].imgH - srcY
                            if (Math.abs(pY) <= sensivity) break;
                        }
                    }
                    if (!this.m_targetTrackClone) {
                        this.m_targetTrackClone = this.createTrackDiv(tracks[src], this.getScreenWidth(), 'track_clone', the_div.dom.childNodes[tracks.length]);
                        this.m_targetTrackClone.dom.style.border = 'green solid 1px';
//                        this.m_targetTrackClone.insertAfter(the_div.dom.childNodes[tracks.length]);
                    }
                    if (dst != src) {
                        if (Math.abs(src - dst) != 1) console.log('src:', src, 'dst:', dst);

                        for (var i = Math.min(src, dst), max = Math.max(src, dst); i < max; i++) this.shiftDownTrack(i);

                        this.m_targetTrackClone.setY(tracks[dst].div.getY() - pY);
                        this.m_XY = xy;
                        // tracks shifting in TrackConfig
                        for (var i = 0; i < tracks.length; i++) {
                            this.m_App.m_Config.TrackConfig[tracks[i].idx].order = i;
                        }
                    } else {
                        this.m_targetTrackClone.setY(tracks[src].div.getY() - delta_y);
                    }
                }
            }
        } else {
             if (e.target.id.indexOf('Ruler') == 1) {
                this.m_vLine.setX(e.getX());
                this.m_vLine.setStyle('display', 'block');
             } else {
                 this.m_vLine.setStyle('display', 'none');
                 if ((e.event.movementY == undefined || (e.event.movementX | e.event.movementY))
                    && Ext.fly(e.getTarget()).hasCls('sv-highlight')) {
                    this.highlightElement(e.getXY());
                }
            }
        }
        e.stopEvent();
    },


    onContextMenu: function(e) {
        e.stopEvent();
        this.showContextMenu(e.getXY());
    },

    // callback allows to customize the menu if you clicked on an object
    showContextMenu: function(xy, callback) {
        if (this.m_Loading) return;
        if (this.m_deferredContext) {
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        var menu = new Ext.menu.Menu();
        var x_pos = xy[0] - Ext.get(this.m_DivId).getX();
        var seq_pos = this.pix2Seq(x_pos - this.m_ScrollPix);
        menu.add(
            {text:'Set New Marker At Position', iconCls:'xsv-markers', scope:this,
             handler:function() {this.pingClick('2-1-0'); this.m_App.newMarkerDlg(this, x_pos); }
        });
        if (this.m_RangeSelectionSet  &&  this.m_RangeSelectionSet.length === 1) {
            menu.add(
                {text:'Set New Marker For Selection',iconCls:'xsv-markers', scope:this,
                 handler:function() {
                    this.pingClick('2-1-1');
                    this.m_App.addMarker(this.m_RangeSelectionSet[0].range);
                    this.removeRangeSelection();
                }
            });
        }
        menu.add('-');
        if (this.m_App.m_Origin) {
            menu.add({text:'Reset Sequence Origin', iconCls:'xsv-origin', scope:this,
                handler:function() {this.m_App.clearOrigin(); this.pingClick('2-1-2');}});
        } else {
            menu.add({text:'Set Sequence Origin At Position', iconCls:'xsv-origin', scope:this,
                handler:function() { this.m_App.setOrigin(this, x_pos); this.pingClick('2-1-2');}});
        }

        if( this.canFlip() ){
            menu.add({
                text:'Flip Sequence Strands', iconCls:'xsv-flip-strands', scope: this,
                enableToggle: true, pressed: this.getFlip(),
                handler:function() { this.flipStrand(); }
            });
        }

        menu.add({iconCls:'xsv-zoom_plus', text:'Zoom In', scope: this,
            handler:function() {this.zoomIn(seq_pos); this.pingClick('2-1-3');}});
        menu.add({iconCls:'xsv-zoom_minus', text:'Zoom Out', scope: this,
            handler:function() {this.zoomOut(seq_pos); this.pingClick('2-1-4');}});
        menu.add({iconCls:'xsv-zoom_seq', text:'Zoom To Sequence', scope: this,
            handler:function() {this.zoomSeq(seq_pos); this.pingClick('2-1-5');}});
        menu.add({iconCls:'xsv-zoom_range', text:'Zoom On Range', scope: this,
            handler:function() {this.zoomRange(); this.pingClick('2-1-6');}});
        if (this.m_App.mf_MultiPanel)
            menu.add({iconCls:'xsv-new_view', text:'Add New Panel on Range', tooltip:'Create New Graphical Panel on Selected Range',
                scope:this, handler:function() {
                    var view = new SeqView.Graphic(this.m_App);
                    this.m_App.registerView(view);
                    var range = this.getTotalSelectedRange();
                    if (range) this.removeRangeSelection();
                    else
                        range = this.m_UrlFrom ? [this.m_UrlFrom, this.m_UrlTo] :
                            [this.m_VisFromSeq + 1, this.m_VisFromSeq + this.m_VisLenSeq];

                    view.startImageLoading(range[0], range[1] - range[0] + 1, {from_ui: true} );
                    this.pingClick('2-1-7');
                }
        });

        menu.add('-');
        // Tools - Blast and Primer Blast
        var noPBlast = (this.m_App.m_ViewParams['acc_type'] !== 'DNA');
        var tools_submenu = new Ext.menu.Menu();
        if (this.m_RangeSelectionSet  &&  this.m_RangeSelectionSet.length > 0  &&  this.m_RangeSelectionSet.length <= 2) {
            tools_submenu.add({iconCls:'xsv-blast', text:'BLAST Search (Selection)', scope: this,
                handler:function() { this.blastSelection(); this.pingClick('2-1-8-1');} });
            tools_submenu.add({iconCls:'xsv-primer', text:'Primer BLAST (Selection)', scope: this, disabled: noPBlast,
                handler:function() { this.primerBlast(); this.pingClick('2-1-8-2');} });
        }
        tools_submenu.add({text:'BLAST Search (Visible Range)', iconCls:'xsv-blast', scope: this, handler:function() {
                              this.m_App.blast([this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1]);
                              this.pingClick('2-1-8-3');
                        } });
        tools_submenu.add({text:'Primer BLAST (Visible Range)', iconCls:'xsv-primer', scope: this, disabled: noPBlast, handler:function() {
                              this.primerBlast(false, [ [this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1] ]);
                              this.pingClick('2-1-8-4');
                        } });
        menu.add({text:'BLAST and Primer Search', tooltip:'Tools', iconCls:'xsv-views_tools', menu:tools_submenu});
        //

        //adding download menu
        var download_submenu = new Ext.menu.Menu();
        download_submenu.add({text:'FASTA (Visible Range)', scope:this,
            handler:function() {this.downloadData(false, "fasta", null); this.pingClick('2-1-9-1');}
        });
        if (this.m_App.m_Markers.length > 0) {
            download_submenu.add({
                text: 'FASTA (All Markers)', scope: this,
                handler: function () { this.downloadAllMarkers("fasta"); this.pingClick('2-1-9-7'); }
            });
        }
        if (this.m_RangeSelectionSet.length > 0) {
            download_submenu.add({text:'FASTA (All Selections)', scope:this,
                handler:function() {this.downloadData(true, "fasta", null); this.pingClick('2-1-9-2');}
            });
        }
        download_submenu.add({text:'GenBank Flat File (Visible Range)', scope:this,
            handler:function() {this.downloadData(false, "flat", null); this.pingClick('2-1-9-3');}
        });
        if (this.m_RangeSelectionSet.length > 0) {
            download_submenu.add({text:'GenBank Flat File (All Selections)', scope:this,
                handler:function() {this.downloadData(true, "flat", null); this.pingClick('2-1-9-4');}
            });
        }

        download_submenu.add(
            {text:'PDF file (Visible Range)', disabled: this.m_App.m_NoPDF, scope: this,
             handler:function() {this.pingClick('2-1-9-5'); this.downloadPDF();}}
        );


        menu.add('-', {text:'Download', iconCls:'xsv-download-static', menu:download_submenu});

        if (this.m_App.m_NoConfDlg !== true && !this.m_App.m_PermConfId) {
            menu.add('-', {iconCls: 'xsv-config', text: 'Configure tracks',
                hidden: this.m_App.countActiveSeqPanels() != 1,
                scope: this, handler:function() { this.m_App.showTracksConfigDlg('2-1-cfg'); }
            });
        }

        if (callback) callback(menu, download_submenu);

        if( menu.items.length > 0 ){
            menu.showAt(xy);
            this.m_ContextMenu = menu;
        }
    },

//////////////////////////////////////////////////////////////////////////
// downloadData:

    downloadData: function(isRangeSelected, format, range) {
        if (isRangeSelected && this.m_RangeSelectionSet.length == 0) {
            Ext.Msg.show({
                title:'Download',
                msg: 'Error: no selection is found',
                buttons: Ext.Msg.OK,
                icon: Ext.Msg.ERROR
           });
           return;
        }
        var ranges = "";
        var first = this.m_App.m_Config.SeqInfo.length, last = 0;
        var flip = this.m_App.m_Flip;
        var strand = flip ? "(-)" : "";

        if (range && range.length > 0) {
            if (!flip)
                ranges = range[0] + "-" + range[1];
            else
                ranges = range[1] + "-" + range[0];
            first = range[0];
            last = range[1];
        }
        else if (!range && isRangeSelected) {
            for (x = 0; x < this.m_RangeSelectionSet.length; x++) {
                range = this.m_RangeSelectionSet[x].range;
                if (x > 0 ) ranges += ",";
                if (!flip)
                    ranges += range[0] + "-" + range[1];
                else
                    ranges += range[1] + "-" + range[0];
                last = Math.max(range[1], last);
                first = Math.min(range[0], first);
            }
        } else {
            last = this.m_VisFromSeq + this.m_VisLenSeq - 1;
            first = this.m_VisFromSeq;
            if (!flip)
                ranges = first + "-" + last;
            else
                ranges = last + "-" + first;
        }
        var url = this.m_App.m_CGIs.SequenceSave + '?id='+ this.m_App.GI + '&format=' + format+ '&ranges='+ ranges
            + '&filename=' + this.m_App.m_Config.SeqInfo.id + '[' + (++first) + '..' + (++last) + ']' + strand + '.' + (format != 'flat' ? 'fa' : format);
        if (this.m_App.m_Key) {
            url += '&key=' + this.m_App.m_Key;
        }
        // Creating form for submitting request to allow browser to handle
        // Content-Disposition header
        if (ranges.length> 0) {
            var form = Ext.DomHelper.append(document.body, {
               tag : 'form',
               method : 'post',
               action : url
           });
           document.body.appendChild(form);
           form.submit();
           document.body.removeChild(form);
        }
    },

    downloadAllMarkers: function (format) {
        if (this.m_App.m_Markers.length == 0) {
            Ext.Msg.show({
                title: 'Download',
                msg: 'Error: no markers found',
                buttons: Ext.Msg.OK,
                icon: Ext.Msg.ERROR
            });
            return;
        }
        var ranges = "";
        var first = this.m_App.m_Config.SeqInfo.length, last = 0;
        var flip = this.m_App.m_Flip;
        var strand = flip ? "(-)" : "";

        for (var m = 0; m != this.m_App.m_Markers.length; m++) {
            var marker = this.m_App.m_Markers[m];
            var range = marker[0];

            if (m > 0) ranges += ",";
            if (!flip)
                ranges += range[0] + "-" + range[1];
            else
                ranges += range[1] + "-" + range[0];
            last = Math.max(range[1], last);
            first = Math.min(range[0], first);
        }
        var url = this.m_App.m_CGIs.SequenceSave + '?id=' + this.m_App.GI + '&format=' + format + '&ranges=' + ranges
            + '&filename=' + this.m_App.m_Config.SeqInfo.id + '[' + (++first) + '..' + (++last) + ']' + strand + '.' + (format != 'flat' ? 'fa' : format);
        if (this.m_App.m_Key) {
            url += '&key=' + this.m_App.m_Key;
        }
        // Creating form for submitting request to allow browser to handle
        // Content-Disposition header
        if (ranges.length > 0) {
            var form = Ext.DomHelper.append(document.body, {
                tag: 'form',
                method: 'post',
                action: url
            });
            document.body.appendChild(form);
            form.submit();
            document.body.removeChild(form);
        }
    },

    downloadPDF: function() {
        var buttHandler = function(butt, event){
            if (butt.itemId != 'btnSave') {
                window.open(this.formPDF.pdf_url + '&inline=true');
                return;
            };
            var form = Ext.DomHelper.append(document.body,
                 { tag: 'form', method: 'post', action: this.formPDF.pdf_url});
            document.body.appendChild(form);
            form.submit();
            document.body.removeChild(form);
        };
        var win = new Ext.Window({
            title: 'Download PDF-file',
            app: this.m_App,
            modal: true,
            layout:'fit',
            minWidth:360,
            width:360,
            height:300,
            plain: true,
            cls: 'SeqViewerApp'
        });
//        if (document.body.clientHeight > win.height) win.renderTo = this.m_App.m_DivId;
        var adjustedRange1=this.m_App.posToLocalDisplay(this.m_VisFromSeq);
        var adjustedRange2=this.m_App.posToLocalDisplay(this.m_VisLenSeq + this.m_VisFromSeq - 1);

        var formPDF = new Ext.FormPanel({
            labelWidth: 1, // label settings here cascade unless overridden
            frame:true,
            bodyStyle:'padding:5px 5px 0',
            pdf_url: '',
            items: [{
                xtype:'fieldset',
                title: 'Enter Sequence Range',
                autoHeight: true,
                defaultType: 'textfield',

                items :[
                    { xtype: 'label', html: 'Possible range formats are '
                        + '10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m<br><br>' },
                    { id: 'sv-pdf_range', width: 290, value: adjustedRange1 + ':' + adjustedRange2,
                        tooltip: 'Range formats are 10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m'
                    },
                    {
                        xtype: 'button', text: 'Create PDF-file',
                        tooltip: 'Range formats are 10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m',
                        width: 300, scope: this, handler:
                            function() {
                                this.m_App.showMessage('', true);
                                this.m_App.showMessage('');
                                this.m_App.m_PDFcomp = Ext.getCmp('sv-PDFcompatibility').checked;
                                this.m_App.m_PDFtitle = Ext.getCmp('sv-PDFtitle').checked;
                                this.m_App.handlePos(Ext.getCmp('sv-pdf_range').getValue(),
                                                     { allow_equal: false, ask_user: true, scope: this,
                                                       success: this.createPDF,
                                                       failure: function(error, options){
                                                           this.m_App.showMessage(error, true);
                                                       }})
                            }
                    }
                ]
            },
            {
            xtype:'checkbox', id: 'sv-PDFtitle', checked: (this.m_App.m_PDFtitle) ? this.m_App.m_PDFtitle : false, boxLabel: 'Add title bar'
            },
            {
            xtype:'checkbox', id: 'sv-PDFcompatibility', checked: (this.m_App.m_PDFcomp) ? this.m_App.m_PDFcomp : true,
                boxLabel: 'Simplified color shading'},{
            xtype:'displayfield', id: 'sv-uplerr' + this.m_App.m_DivId, value: '', hidden: true,
                style: {color:'red', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}},{
            xtype: 'displayfield', id: 'sv-uplmsg' + this.m_App.m_DivId, value: 'Applying simplified color shading provides wider compatibility with 3rd party PDF-editors',
                style: {color:'grey', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}}
            ],

            buttons: [{
                text: 'View', itemId: 'btnView', scope: this, disabled: true, handler: buttHandler },{
                text: 'Save', itemId: 'btnSave', scope: this, disabled: true, handler: buttHandler },{
                text: 'Cancel', scope: this, handler: function(){ win.close(); delete this.formPDF; }
            }]
        });
        win.add(formPDF);
        this.formPDF = formPDF;
        win.show();
    },

    setPreferences: function() {
        var app = this.m_App;
        var optWindow = new Ext.Window({
            title: 'Preferences',
            app: app,
            layout: 'form',
            modal: true,
            buttons: [{text: 'Update',
                handler: function() {
                    var preps = {};
                    Ext.each(optWindow.items.items, function(p){preps[p.getName()] = p.getValue()});
                    Ext.apply(app.m_Config.Options, preps);
                    localStorage.setItem('NCBI/SV/Preferences', window.JSON.stringify(preps));
                    optWindow.close();
                    app.reloadAllViews();
                }},
                {text: 'Cancel', handler: function() {optWindow.close();}}]
        });
//        if (document.body.clientHeight > optWindow.height) optWindow.renderTo = this.m_App.m_DivId;

        var fields = [['color','Coloration'], [ 'label','Label Placement'], ['decor','Shapes'], ['spacing','Spacing']];
        var options = app.m_Config.Options;
        Ext.each(fields, function(p) {
            optWindow.add({
                xtype: 'combo', triggerAction: 'all', fieldLabel: p[1], mode: 'local', name: 'curr_' + p[0],
                store: options.controls[p[0]], allowBlank: false, editable: false, value: options['curr_' + p[0]]
            });
        });
        optWindow.show();
    },

    createPDF: function(pos_range, options){
        var width = this.getScreenWidth() * (pos_range[1] - pos_range[0] + 1) / this.m_VisLenSeq;
        if (width > 8000){
            this.m_App.showMessage('Image is too wide (projected  width is ' + width + ' which is more than 8000 permitted)', true);
            return;
        }
        this.formPDF.getEl().mask('Creating PDF ...');
        this.formPDF.pdf_url = '';
        this.formPDF.down('#btnView').disable();
        this.formPDF.down('#btnSave').disable();
        // Instead of parsing/unparsing markers and making coords zero-based we instruct getGraphicParams
        // to return us stuff prepared for CGI use.
        var params = this.getGraphicParams(pos_range[0], pos_range[1] - pos_range[0] + 1, width, true);
        params.simplified = this.m_App.m_PDFcomp;
        params.target = 'pdf';
        params.tracks = params.tracks[0];
        params.pdftitle = this.m_App.m_PDFtitle;
        this.formPDF.filename = this.m_App.m_Config.SeqInfo.id + '[' + (params.from + 1) + '..' + (params.from + params.len) + ']'
        this.m_App.AjaxRequest({url: this.m_App.m_CGIs.Graphic, context: this, data: params,
            success: this.checkJobStatusPDF,
            error: this.checkJobStatusPDF
        });
    },

    checkJobStatusPDF: function(data, text, res) {
        var from_cgi = SeqView.decode(data);
        if (from_cgi.job_status) {
            if (from_cgi.job_status == 'failed' || from_cgi.job_status == 'canceled') {
                this.formPDF.getEl().unmask();
                this.m_App.showMessage(from_cgi.error_message, true)
            } else {
                if (from_cgi.progress_message && from_cgi.progress_message.length > 0) {
                    var progress_text = from_cgi.progress_message.replace(/\&quot;/gi, "\"");
                    var current_task = Ext.decode(progress_text).current_task;
                    if (current_task)
                        this.m_App.showMessage(current_task);
                }
                Ext.defer(SeqView.App.simpleAjaxRequest, 500, this, [{
                    url: this.m_App.m_CGIs.Graphic + '?job_key=' + from_cgi.job_id,
                    context: this,
                    success: this.checkJobStatusPDF,
                    error: this.checkJobStatusPDF }]);
            }
       } else {
            if (from_cgi.pdf_url){
                // If the pdf_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
                // This is a way to provide reliable URL resolution for embedding. SV-1760
                if (from_cgi.pdf_url && from_cgi.pdf_url.charAt(0) == '?') {
                    from_cgi.pdf_url = this.m_App.m_CGIs.NetCache + from_cgi.pdf_url;
                }
                this.formPDF.pdf_url = from_cgi.pdf_url + '&filename=' + this.formPDF.filename + '.pdf';
                this.formPDF.down('#btnView').enable();
                this.formPDF.down('#btnSave').enable();

                var scope = this;
                var url = this.formPDF.pdf_url + '&inline=true';
                NCBIGBUtils.makeTinyURL(url, function(res) {
                    if (res.id) url = res.id;
                    scope.formPDF.getEl().unmask();
                    scope.m_App.showMessage('<a href=' + url +' target=\"_blank\">' + url + '</a>');
                    scope.formPDF.updateLayout();
                });
            } else {
                this.formPDF.getEl().unmask();
                this.m_App.showMessage('Request failed', true);
            }

        }
    },


    onDblClick: function(e) {
        if (!Ext.fly(e.getTarget()).hasCls('sv-dblclick'))
            return;

        if (this.m_Loading)
            return;

        var the_elem_xy = Ext.get(this.m_DivId).getXY();
        var xx = e.getXY()[0] - the_elem_xy[0] - this.m_ScrollPix;

        var pix = Math.abs(xx);//config['scroll_pix']) + e.getXY()[0] + the_elem_xy[0];
        var seq_pos = this.pix2Seq(pix);

        this.zoomIn(seq_pos);
    },

    createMarkerElem: function(marker) {
        var elem = Ext.get(this.m_DivId);
        var create_params = marker.getCreateParams(false, this.m_Idx);
        return create_params.template.append(elem, create_params.options, true);
    },

    changeSelectedSig: function(area, ctrl_key) {
        // don't do anything for "graphics" and "histograms" and for fake areas
        if (area.type === undefined || (area.type & SeqView.AreaFlags.NoSelection) !== 0)
            return;

        var cur_sig = this.m_SelectedSig ? this.m_SelectedSig.split(';') : [];
        var new_sig = area.signature;
        if ( ctrl_key ) {
            var s_idx = cur_sig.indexOf(new_sig);
            if (s_idx != -1) {
                cur_sig.splice(s_idx,1);
            } else {
                cur_sig.push(new_sig);
            }
            new_sig = cur_sig.join(';');
        } else {
            this.removeSelectionsWithNoPinnedToolTips();
        }
        if (this.setSelectedSig(new_sig)) {
            this.pingClick('6-0');
            this.m_App.fireEvent('feature_clicked', this, area);
        }

    },

    setSelectedSig: function(new_sig) {
        if (!this.m_SelectedSig && !new_sig) return false;

        this.m_SelectedSig = (this.m_SelectedSig === new_sig) ? "" : new_sig;
        this.refresh();
        return true;
    },

    zoomRange: function() {
        var range = this.getTotalSelectedRange();
        if (range) {
            this.removeRangeSelection();
            this.startImageLoading(range[0], range[1] - range[0] + 1, {from_ui: true} );
        } else  this.gotoPositionDlg(true);
    },

    zoomIn: function(seq_pos) {
        if (seq_pos == undefined && this.getMinLength() >= this.m_VisLenSeq) return;
        var new_len = Math.floor(this.m_VisLenSeq / 2);
        var center = seq_pos? seq_pos : this.m_VisFromSeq + new_len;
        var new_from = Math.floor(Math.max(0, center - new_len/2));
        if (new_from + new_len > this.m_App.m_SeqLength) new_len = this.m_App.m_SeqLength - new_from;

        this.startImageLoading(new_from, new_len, {from_ui: true});
    },

    zoomSeq: function(center_seq_pos){

        if( !center_seq_pos ){
            center_seq_pos = Math.floor(this.m_VisFromSeq + this.m_VisLenSeq / 2);
        }

        var new_len  = this.getMinLength();;
        var new_from = Math.floor(center_seq_pos - new_len / 2);
        if( new_from < 0 ){
            new_from = 0;
        }
        var app = this.m_App;
        if( new_from + new_len > app.m_SeqLength ){
            new_len = app.m_SeqLength - new_from;
        }
        var tracks_changed = app.showTracks({key: 'sequence_track'}, true);
        this.startImageLoading( new_from, new_len, {from_ui: true} );
        if (tracks_changed) {
            app.fireEvent('configuration_changed', app);
        }
    },

    zoomOut: function(seq_pos) {
        if (seq_pos == undefined && this.m_VisLenSeq >= this.m_LenSeq) return;
        var new_len  = this.m_VisLenSeq * 2;
        var center = seq_pos? seq_pos : this.m_VisFromSeq +  this.m_VisLenSeq / 2;
        var new_from = Math.max(0, center - new_len/2);
        if (new_from + new_len > this.m_App.m_SeqLength) new_len = this.m_App.m_SeqLength - new_from;
        this.startImageLoading( new_from, new_len, {from_ui: true} );
    },

    parseAndGotoPosition: function(text, range_only) {
        this.m_App.parseAndGotoPosition(text, {range_only: range_only, view: this});
    },

    gotoAndSearch: function(term) {
        var combo = this.m_View.down('#gotoBox');
        if (term == null || term.length == 0) {
            if (combo) {
                term = combo.getValue();
                combo.setValue('');
            }
       }
        if (!term) return;

        this.m_App.addSearchPattern(term);
        if (combo) combo.setStore(this.m_App.searchPatternData);
        this.m_App.saveSearchPatternData();

        this.m_App.gotoAndSearch(term, {view: this});
    },

    gotoPosRange: function(pos_range, center, options ){
        var new_from = pos_range[0];
        var len = this.m_VisLenSeq;
        if (pos_range.length === 2 && pos_range[1]) { // range is specified
            len =  pos_range[1] - new_from + 1;
        } else { // single position
            // according to SV-620 we want to keep the length of the current view in case if
            // it is less then 80% (?? it may change) of the sequence length
            if (len >= this.m_App.m_SeqLength * 0.8 ) {
                var new_center = new_from;
                if (new_center - 500000 < 0 && new_center + 500000 > this.m_App.m_SeqLength) {
                    if (new_center - 500 < 0 && new_center + 500 > this.m_App.m_SeqLength) {
                        var screen_width = this.getScreenWidth();
                        new_from = new_center - screen_width/8/2;
                        len = screen_width/8;
                    } else {
                        new_from = new_center - 500;
                        len = 1000;
                    }
                } else {
                    new_from = new_center - 500000;
                    len = 1000000;
                }
            } else if (center) new_from -= len / 2;

            if (new_from < 0) new_from = 0;
            if (new_from + len > this.m_App.m_SeqLength) {
                new_from -= new_from + len - this.m_App.m_SeqLength;
                if (new_from < 0) {
                    new_form = 0;
                    len = this.m_App.m_SeqLength;
                }
            }
        } // end single position
        this.startImageLoading(new_from, len, options);
    },

    gotoPositionDlg: function(range_only) {
        Ext.MessageBox.prompt(range_only ? "Zoom to range" : "Go to position/range",
          'Enter sequence ' + (range_only ? '' : 'position or') +
              ' range<br />(possible range formats are 10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m):',
          function(btn, text) {
              if (btn == 'ok') this.parseAndGotoPosition(text, range_only);
          }, this, false);
    },

    primerBlast: function(whole, range_set) {
        //check that range is not exceeding 50Kb

        var seqinfo = this.m_App.m_ViewParams;
        if (seqinfo["acc_type"] !== "DNA") return;
        range_set = range_set || this.getSelectedRangeSet();
        if (!whole  && (!range_set || range_set.length < 1)) {
            //open message window
            Ext.MessageBox.show({
               title: 'Primer BLAST (Selection)',
               msg: 'Cannot run BLAST - select a range first',
               buttons: Ext.MessageBox.OK//,
           });
           return;
        } else if (!whole  && (range_set && range_set.length > 2)) {
            //open message window
            Ext.MessageBox.show({
               title: 'Primer BLAST (Selection)',
               msg: 'Cannot run Primer BLAST on more than 2 ranges',
               buttons: Ext.MessageBox.OK//,
           });
           return;
        } else if (range_set && range_set.length == 1){
            var range = range_set[0];
            if(range[1]-range[0]>50000) {
                Ext.MessageBox.show({
               title: 'Primer BLAST',
               msg: 'Cannot run Primer BLAST on a range longer than 50K',
               buttons: Ext.MessageBox.OK//,
               });
               return;
            }
        } else if (range_set && range_set.length == 2){
            var range1 = range_set[0];
            var range2 = range_set[1];
            if((range2[1]-range1[0])>50000) {
                Ext.MessageBox.show({
               title: 'Primer BLAST',
               msg: 'Cannot run Primer BLAST on ranges longer than 50K',
               buttons: Ext.MessageBox.OK//,
               });
               return;
            }
        }
        if (whole)
            range_set = [ [0, this.m_App.m_SeqLength-1] ];
        this.m_App.primerBlast(whole, range_set);
    },

    blastSelection: function() {
        if (!this.m_RangeSelectionSet || this.m_RangeSelectionSet.length < 1)
        {
            Ext.MessageBox.show({
               title: 'BLAST Search (Selection)',
               msg: 'Cannot run BLAST - select a range first',
               buttons: Ext.MessageBox.OK
           });
            return;
        }
        var range = this.getTotalSelectedRange();
        this.m_App.blast(range);
    },

    openFullView: function(mode) {
        if (mode !== 'full') mode = 'portal';
        this.m_App.getLinkToThisPageURL(mode);
    },

    showAssmInfo: function() {
        var sts = ['Not set in the URL (and could not be guessed)',
                'Explicitly set from the URL',
                'Automatically determined by the provided molecule ID',
                'Detection failed due to timeout',
                'Conflict between the URL parameters: molecule ID and Assembly context',
                'Automatically determined by the provided molecule ID (using cache)',
                'Obtained from the molecule record'];
        var seqInfo = this.m_App.m_Config.SeqInfo;
        var msg = 'Name: <a href="' + SeqView.webNCBI + 'assembly/' + seqInfo.assm_context
            + '/" target=\"_blank\">' + seqInfo.assm_info.name
            + '</a><br>TaxID: <a href="' + SeqView.webNCBI + 'Taxonomy/Browser/wwwtax.cgi?id=' +
            + seqInfo.assm_info.tax_id + '" target=\"_blank\">' + seqInfo.assm_info.tax_id
            + '</a><br> Type: ' + (seqInfo.assm_info.type ? 'GenBank' : 'RefSeq')
            + '<br>' + sts[seqInfo.assm_context_status];
        if (seqInfo.warnTT) msg += '<br>WARNING!<br>' + seqInfo.warnTT;
        Ext.MessageBox.show({
            icon: seqInfo.warnTT ? Ext.MessageBox.WARNING : Ext.MessageBox.INFO,
            title: 'Genomic Assembly Context',
            msg: msg,
            buttons: Ext.MessageBox.OK
        });
        return;
    },

    setGeneMode: function(mode, config) {
        var bttn = this.m_View.down('#geneModeButton');
        if (!bttn) return;
        config = config || this.m_App.m_Config;
        if (this.m_showAllGenes != mode) config.changeGeneMode(this.m_showAllGenes = mode);
        bttn.setPressed(mode);
        bttn.setTooltip('Switch ' + (bttn.pressed ? 'OFF' : 'ON') + ' mode "show All" for Gene tracks');
        return bttn;        
    },

    hideGeneModeButton: function(hide) {
        var bttn = this.setGeneMode(hide ? false : this.m_showAllGenes || false);
        if (bttn) bttn[hide ? 'hide' : 'show']();    
    }
});  // End Of SeqView.Graphic
/*  $Id: tm_panels.js 38325 2017-04-25 22:08:59Z borodine $
*/
Ext.namespace('SeqView.TM');

SeqView.TM.updateTrackDetails = function(panel, track){
    if (!track.choice_list && !track.check_boxes) {
        panel.add({xtype:'displayfield', hideLabel: true, value: "<i>No settings available</i>"});
        return;
    }
    var fPanel = new Ext.Panel({
        border: false,
        bodyStyle: { background: 'inherit' }, //'#F0F0F0'},
        layout: 'form',
        labelWidth: 130
    });

//    var genesRenderMode = '';
    if (track.choice_list) {
/*        if (track.category.name == 'Genes') {
            genesRenderMode = track.choice_list[0].curr_value;
        }*/
        Ext.each(track.choice_list, function(lst) {
            var data = [];
            var cur_text = '';
            Ext.each(lst.values, function(v) {
                data.push({text: v.display_name, value: v.name, help: v.help});
                if (v.name == lst.curr_value) cur_text = v.display_name;
            });
            var store = Ext.create('Ext.data.Store', {
                model: Ext.define(null, { extend: 'Ext.data.Model',
                   fields: ['text', 'value', 'help']}),
                data: data });
            var cb = Ext.create('Ext.form.field.ComboBox', {
                allowBlank: false,
                choiceLst: lst,
                disabled: !track.shown,                
                queryMode: 'local',
                store: store,
                editable: false,
                fieldLabel: lst.display_name,
//                tooltip: {text: lst.help},
/*                listConfig: {
                    tpl: [
                        '<table width="100%"><tpl for=".">',
                            '<tr data-qtip="{help}">',
                                '<td class="x-boundlist-item">{text}</td>',
                            '</tr>',
                        '</tpl></table>'
                    ]
                },*/

/*                listConfig: {
                    getInnerTpl: function() {
                        return '<div data-qtip="{help}">{text}</div>';
                    }
                },*/
                listeners: {
                    beforeselect: function(f, r) {
/*                        if (genesRenderMode) {
                            panel.changedGenesMode |= genesRenderMode != r.data.value;
                        }*/
                        f.choiceLst.curr_value = r.data.value;
                    }
                },
                value: cur_text,
                width: 250
            });
            fPanel.add(cb);
        });
    }
    if (track.check_boxes) {
        for(var i = 0; i < track.check_boxes.length; ++i) {
            var cbox = track.check_boxes[i];
            var item = new Ext.form.Checkbox({
                boxLabel: cbox.display_name,
                checked: cbox.value,
                disabled: !track.shown,
                tooltip: cbox.help,
                rcbox: cbox,
                listeners: { change: function(cb, newVal) { cb.rcbox.value = newVal; } }
            });
            if (i == 0) {
                item.fieldLabel = 'Other Settings';
            } else {
                item.fieldLabel = '&nbsp;';
                item.labelSeparator = '';
            }
            fPanel.add(item);
        }
    }
    if (track.text_boxes)
        for (var i = 0; i < track.text_boxes.length; i++) {
            var tbox = track.text_boxes[i];
            if (tbox.name == 'HistThreshold') {
                var item = new Ext.form.NumberField({
                    fieldLabel: tbox.display_name,
                    allowNegative: false,
                    disabled: !track.shown,
                    emptyText: 'None',
                    value: tbox.value || 0,
                    rtbox: track.text_boxes[i],
                    tooltip: tbox.help,
                    listeners: {
                        change: function(tb, nval, oval) {
                            if (nval == 0) tb.setValue('undefined');
                            tb.rtbox.value = nval || 0;
                        },
                        render: function(tb) {
                            if (tb.value == 0) tb.setValue('undefined');
/*                            Ext.QuickTips.register({
                                target: tb.getEl(),
                                text: tb.tooltip
                            });*/
                        }
                    }
                });
            }
            if (item) fPanel.add(item);
        }
    panel.add(fPanel);
    panel.updateLayout();
};

Ext.define('SeqView.TM.TrackDetail', {
    extend: 'Ext.Panel',
    border: false,
    bodyStyle: {background: 'inherit'}, //'#F0F0F0'},
//    itemId: 'track_details',
    overflowY: 'auto',
 
    title: 'Track Settings',
    updateDetail: function(rdata) {
        var panel = this;
        this.items.each(function(item) {
            panel.remove(item);
        });
        panel.setTitle("Track Settings: " + rdata.display_name || rdata.name);
        if (rdata.help) {
            var ref  = !rdata.legend_text ? ''
                : '<br/><a href=\"' + SeqView.getHelpURL() + 'legends/#' + rdata.legend_text
                + '\" target=\"_blank\" style=\"color:blue\">Track legend</a>';
            panel.add({ xtype:'fieldset', padding: 7, html: rdata.help + ref });
        }
        SeqView.TM.updateTrackDetails(panel, rdata);
        panel.ownerCt.updateLayout();
    }
});

SeqView.TM.collectAllTracks = function(categories, activeOnly) {
    var collected_tracks = [];
    Ext.each(categories, function(cat) {
        Ext.each(cat.subcategories, function(subcat) {
            Ext.each(subcat.tracks, function(track) {
                if (track.shown || !activeOnly ) {
                    if (cat.name && cat.name != 'other')
                        track.cat_name = cat.name;
                    else
                        delete track.cat_name;
                    if (subcat.name && subcat.name != 'other')
                        track.subcat_name = subcat.name;
                    else
                        delete track.subcat_name;
                    collected_tracks.push(track);
                }
            });
        });
    });
    collected_tracks.sort(function(t1,t2) {
        return t1.order - t2.order;
    });
    return collected_tracks;
}

SeqView.TM.collectActiveTracks = function(categories) {
    return SeqView.TM.collectAllTracks(categories, true);
}


Ext.define('SeqView.TM.TracksPanel', {
    extend: 'Ext.Panel',
    initComponent: function() {
        if (this.activeTracks) {
            this.gridPanel = new SeqView.TM.GridBase({ tracks: this.tracks, sortable: false,
                 multiSelect: true,
                 viewConfig: { plugins: { ptype: 'gridviewdragdrop', containerScroll: true },
                     listeners: { drop: function() {
                         var tracks = this.grid.tracks;
                         var order = 0;
                         this.grid.getStore().data.each(function(r) {
                              Ext.each(tracks, function(t) {
                                  if (t.uuid == r.data.uuid) {
                                       t.order = order++;
                                       return false;
                                  }
                              });
                         });
                     }}
                 }});
        } else 
            this.gridPanel = (this.searchTracks)
                ? SeqView.TM.createSearchTracksTab(null, this.tracks_config)
                : new SeqView.TM.GridBase({tracks: this.subcategories ? this.subcategories[0].tracks : this.tracks,
                                                        subcategories: this.subcategories});

        if (!this.searchTracks) {
            this.detailPanel = new SeqView.TM.TrackDetail({
                region: 'south',
                split: true,
                height: 250,
                minSize: 100,
                collapsible: true
            });
            var sm = this.gridPanel.getSelectionModel();
//            if (this.uploadedData) sm.setSelectionMode('MULTI');
            sm.on('selectionchange', this.onRowSelect, this);
        }
        var items = [];
        if (this.activeTracks) {
            var hintStr = {
                region: 'north',
                border: false,
                bodyStyle: { padding: '5px 5px 5px 5px', 'font-size': '14px' },
                html: 'Click on track to display settings. To re-order tracks, drag and drop track names.'
            };
            items.push(hintStr, this.gridPanel, this.detailPanel);
        } else if (this.searchTracks) {
            var sbar = new Ext.Panel({
                region: 'north',
                border: false,
                tbar: ['-',' Search: ', ' ', new SeqView.SearchField({store: this.gridPanel.store, width: 220})]
            }); 

            items.push(sbar, this.gridPanel);
        }
        else items.push(this.gridPanel, this.detailPanel);
        Ext.apply(this, {
            border: false,
            layout: 'border',
            items: items
        });
        this.callParent(arguments);

        this.gridPanel.on('viewready', function(g) {
            g.getSelectionModel().select(0);
        });
    },

    onRowSelect: function(sm, srow) {
        if (!srow.length) return;
        var rdata = srow[0].getData()
        this.detailPanel.updateDetail(rdata);
        var group_button = Ext.getCmp('sv-config_group_button_id');
        var group = rdata.setting_group;
        if (group) {
            group_button.show();
            group_button.track = rdata;
            var display_name = group;
            var display_limit = 10;
            if (group.length > display_limit) {
                display_name = group.substr(0, display_limit) + "...";
            }
            group_button.setText("Configure Group '" + display_name + "'");
        } else {
            group_button.hide();
        }
    }
});

Ext.define('SeqView.TM.MainPanel', {
    extend: 'Ext.Panel',
    initComponent : function(){
        Ext.apply(this, {
            title: 'Tracks',
            layout: 'fit',
            itemId: 'trackpanel',
            dirty: false,
            items: this.createPanel(),
            listeners: {
                show: function(p) {
                    p.updateLayout();
                }
            }
        });
        this.callParent(arguments);
    },

    reload: function(categories) {
        this.removeAll();
        this.categories = categories;
        this.dirty = false;
        var panel = this.createPanel();
        this.add(panel);
        if (this.isVisible())
            this.updateLayout();
    },

    createPanel: function() {
        var tracks_config = this.sviewApp.m_Config.TrackConfig;
        var categories = this.categories;
        var tabs = [];
        this.activeTracksTab = new SeqView.TM.TracksPanel({
            title: 'Active Tracks',
            iconCls: 'xsv-seq-logo',
            activeTracks: true,
            tracks_config: tracks_config,
            tracks: SeqView.TM.collectActiveTracks(categories),
            listeners: {
                deactivate: function() { this.gridPanel.syncModifiedTracks(this.tracks_config); }
            }
        });
        tabs.push(this.activeTracksTab);
        this.searchTracksTab = new SeqView.TM.TracksPanel({
            title: 'Search Tracks',
            iconCls: 'xsv-search-button',
            searchTracks: true,
            tracks_config: tracks_config,
            searchStr:'*',
            listeners: {
                deactivate: function(p) { this.gridPanel.syncModifiedTracks(); }
            }
        });
        tabs.push(this.searchTracksTab);
        for (var i = 0; i < categories.length; i++) {
            var cat = categories[i];
            var noUploads = true;
            Ext.each(cat.subcategories, function() {
                Ext.each(this.tracks, function(){ return noUploads = this.id.charAt(0) != 'U'; });
                return noUploads;
            });
            var catTab = new SeqView.TM.TracksPanel({
                title: cat.display_name,
                layout: 'fit',
                category: cat,
                subcategories: cat.subcategories,
                tracks_config: tracks_config,
                uploadedData: !noUploads,
                listeners: {
                    activate: function(p) {
                        if (this.uploadedData) Ext.getCmp('sv-delete_button_id').enable(); 
                    },
                    deactivate: function(p) {
                        if (this.uploadedData) Ext.getCmp('sv-delete_button_id').disable(); 
                        this.gridPanel.syncModifiedTracks(this.tracks_config);
                    }
                }
            });
            tabs.push(catTab);
        }
        var panel = new Ext.tab.Panel({
            activeTab: 0,
            tabPosition: 'left',
            tabRotation: 0,
            plain:true,
            tabWidth: 158,
            items: tabs
        });
        return panel;
    },
    
    getActiveGridPanel: function() {
        return this.items.items[0].getActiveTab().gridPanel;  
    },

    setDirty: function(flag) {
        this.dirty = (flag === true);
    }
});

Ext.define('SeqView.TM.UploadPanel', {
    extend: 'Ext.Panel',
    initComponent : function(){
        Ext.apply(this, {
            title: 'Custom Data',
            layout: 'border',
            itemId: 'uploadpanel',
            infoMsg: 'No informations or details',
            items:[{
                region: 'west', layout:'fit', width: 168, title:'Data Source',
                items: [{
                    xtype: 'treepanel', itemId: 'treepanel', rootVisible: false, enableDD: false, lines: false,  autoScroll: false,
                    listeners: { itemClick: function(node, rec) {
                        this.uploadButt.disable();
                        this.updateDiffPart(rec.data.type, rec.data.text);
                    }, scope: this},
                      root: {
                        children: [
                            {text:'BLAST Results', type:'rid', iconCls:'xsv-ext_data', leaf:true},
                            {text:'Data File', type:'file', iconCls:'xsv-ext_data', leaf:true},
                            {text:'BAM Files', type:'BAMfiles', iconCls:'xsv-ext_data', leaf:true},
                            {text:'Alignment MUSCLE/FASTA', type:'FASTAfile', iconCls:'xsv-ext_data', leaf:true},
                            {text:'URL', type:'url', iconCls:'xsv-ext_data', leaf:true},
                            {text:'Text', type:'text', iconCls:'xsv-ext_data', leaf:true}
                    ]}
                }]
            },{
                frame:true,
                xtype:'form',
                region: 'center',
                hidden: false,
                itemId:'loaddatapanel',
                fileUpload: true,
                items:[                        
                    {xtype:'panel', layout: 'form', border: false, bodyStyle: {background: 'inherit'},
                    items:[
                        {xtype:'textfield', fieldLabel: 'Track Name', name:'track_name'},
                        {xtype:'hidden', name:'assm_acc', value: this.cfgPanel.sviewApp.m_AssmContext || ""}]},
                    {xtype:'displayfield', id: 'sv-uplmsg' + this.cfgPanel.sviewApp.m_Idx, value: '', width:'95%',
                        style: {color:'grey', 'white-space': 'nowrap'}},
                    {xtype:'displayfield', id: 'sv-uplerr' + this.cfgPanel.sviewApp.m_Idx, value: '',
                        style: {color:'red', 'white-space': 'nowrap'}},
                    {xtype: 'button', hidden: true, scope: this, handler: this.showErrorDetails,
                        id: 'sv-err_details_button_id'}
                ]
            }],
            buttons: [{
                text:'Upload', scope: this, itemId: 'uploadButton',
                disabled: true,
                handler: function(b, e) {
                    e.stopEvent();
                    var form = this.down('#loaddatapanel').getForm();
                    if (!form.isValid()) return;
                    if (form.isDirty() || this.dzfiles) {
                        this.submitData(form);
                    }
                    else {
                        this.cfgPanel.sviewApp.showMessage("Inputs are empty or invalid!", true);
                    }
                    this.cfgPanel.sviewApp.pingClick('9-5-' + this.type);
                }
            }]
        });
        this.callParent(arguments);
        this.down('#treepanel').getSelectionModel().select(0);
        this.updateDiffPart('rid', 'BLAST Results');
        this.uploadButt = this.down('#uploadButton');
    },
    showErrorDetails: function() {
        if (this.infoMsg == '') return;
        var mbw = new Ext.Window({
            title: 'Error details',
//            app: this.cfgPanel.sviewApp,
            modal: true, layout:'fit',
            width: 400, height: 300});
        mbw.add(new Ext.FormPanel({ labelWidth: 1, autoScroll: true,
            items:[
               {xtype: 'displayfield', value: this.infoMsg, textalign: 'left'}],
        	   buttons:[{text: 'OK', handler: function() {mbw.close();}}]})
        );
        mbw.show();
    },
    Sec2Time: function(sec) {
        var hours   = Math.floor(sec / 3600);
        var minutes = Math.floor((sec - (hours * 3600)) / 60);
        var seconds = sec - (hours * 3600) - (minutes * 60);
        var time = '';
        if (hours > 0) 
            time = hours + ' hour' + ((hours > 1) ? 's ' : ' ');
        if (minutes > 0)
            time += minutes + ' minute' + ((minutes > 1) ? 's ' : ' ');
        if (seconds > 0)
            time += seconds + ' second' + ((seconds > 1) ? 's' : '');
        return time;
    },
    updateMessage: function() {
        this.currTime++;
        this.totalTime++;

        var msg = '<br>';
        msg += 'Total time: ' + this.Sec2Time(this.totalTime) + '<br>';
        for (var i = 0; i < this.tasks.length; i++)
            msg += this.tasks[i].task + ": " + this.tasks[i].time + " seconds<br>";
        msg += this.currTask + ": " + this.Sec2Time(this.currTime);
        if (this.percentage) {
            var time = this.Sec2Time(Math.round(this.currTime * (100/this.percentage - 1)));
            if (time.length > 0) 
                msg += " (" + time + " remaining)";
        }
        msg += '<br>';
        if (this.percentage) msg += 'Percentage: ' + this.percentage + '%<br>';
        this.cfgPanel.sviewApp.showMessage(msg);
        Ext.defer(this.updateMessageWrap, 1000, this, []);
    },

    updateMessageWrap: function() {
        if (this.uploaderUUD) this.updateMessage("");
    },

    reset: function(panel) {
        if (this.dzfiles) {
            for (var f in this.dzfiles) {
                this.dzfiles[f].dzfname.update('');
                this.dzfiles[f].inputFld.value = '';
            }
            delete this.dzfiles;
        } else 
            panel.getForm().setValues({rid:'', data: '', dataURL: ''});
    },

    cleanupUpload: function(msg) {
        delete this.uploaderUUD;

        var button = Ext.getCmp('sv-configure_button_id')
        if (button) button.enable();
        if (msg) this.cfgPanel.sviewApp.showMessage(msg, true);
    },
    
    submitData: function(form, callback, scope) {
        this.currTime = this.totalTime = this.percentage = 0;
        this.currTask = 'uploading file';
        this.tasks = [];

        Ext.getCmp('sv-err_details_button_id').hide();
        Ext.getCmp('sv-configure_button_id').disable();
        this.uploadButt.disable();
        this.cfgPanel.sviewApp.showMessage('');
        this.cfgPanel.sviewApp.showMessage('', true);

        var config = {assm_acc: this.cfgPanel.sviewApp.m_AssmContext};
        var fval = form.getValues();
        for (var v in fval) {
            var val = fval[v];
            if (typeof val == 'string' && val.length == 0) continue;
            config[v] = val;
        }
        var find_comp = config.find_comp || false;
        delete config.find_comp;
        if (config.rid) config.blast = {rid: config.rid, link_related_hits: find_comp};
        else config.check_cs = true;
        var uPanel = this;
        if (callback) {
            this.exitCallback = {callback: callback, scope: scope};
        }
        uPanel.consError = console.error;
        console.error = function() {}
        var finalize = function(msg) {
            console.error = uPanel.consError || console.error;
            delete uPanel.consError;
            uPanel.cleanupUpload(msg);
            if (!msg) {
                var msg = 'Data uploaded';
                if (uPanel.cfgPanel.sviewApp.m_AssmContext) {
                    msg += " on assembly " + uPanel.cfgPanel.sviewApp.m_AssmContext;
                }
                uPanel.cfgPanel.sviewApp.showMessage(msg);
                var dTrk = Ext.each(uPanel.infoUpTracks, function() { return (this.track_type == 'non-displayable_track') });
                if (typeof dTrk == 'undefined') uPanel.cfgPanel.sviewApp.showMessage('There is no data that can be displayed on the sequence', true);
            }
            uPanel.reset(uPanel.down('#loaddatapanel'));
            if (uPanel.exitCallback) {
                var tmpobj = uPanel.exitCallback;
                delete uPanel.exitCallback;
                tmpobj.callback.call(tmpobj.scope);
            }
        }

        if (this.dzfiles) Ext.apply(config, this.dzfiles);
        try {
            this.uploaderUUD = new UUD.FileUploader(config);
            var promise = this.uploaderUUD.getPromise();
            promise.fail(function(){
                finalize('Failed to upload data: ' + this.getErrors());
            });
            promise.done(function(tlist, dkey) {
                uPanel.infoUpTracks = this.getTracks();
                if (uPanel.infoUpTracks) {
                    Ext.getCmp('sv-load_default_button_id').disable();
                    Ext.each(uPanel.infoUpTracks, function() { uPanel.cfgPanel.sviewApp.addUploadedTrackID(this);});
                    uPanel.ownerCt.down('#trackpanel').setDirty(true);
                }
                var errMsg = this.getErrors();
                if (errMsg.length) {
                    uPanel.infoMsg = '';
                    Ext.each(errMsg, function(msg, idx){
                        this.infoMsg +='# ' + (idx + 1) + '. ' + msg + '<br>';
                    }, uPanel);
                    
                    var bttn = Ext.getCmp('sv-err_details_button_id');
                    var bttxt ='Data parsing error details (' + errMsg.length + ')';
                    bttn.setText(bttxt);
                    bttn.show();
                }
                finalize();
            });
            promise.progress(function(progress) {
                if (progress.percentage) uPanel.percentage = progress.percentage;
                var task = progress.current_task;
                if (task == "" || task == 'pending' || task == uPanel.currTask) return;

                uPanel.tasks.push({task: uPanel.currTask, time: uPanel.currTime});
                uPanel.currTask = task;
                uPanel.currTime = 0;
            });
            this.uploaderUUD.upload();
        } catch(e) { finalize('Unable to upload data: ' + e.message); } 
    
        this.updateMessage();
    },

    updateDiffPart: function(type, title) {
        var uPanel = this;
        this.type = type;
        var fileUploadDisabled = Ext.isIE && !(Ext.isIE11 || Ext.isIE10);
        var manageButt = function(str) {
            uPanel.uploadButt[str.length ? 'enable' : 'disable']();
            return true;
        };
        var processFName = function(ffld, file) {
            file.inputFld = document.getElementById(ffld.getInputId());
            file.dzfname = uPanel.down('#dz_' + ffld.itemId);
            uPanel.dzfiles = uPanel.dzfiles || {};
            uPanel.dzfiles[ffld.itemId] = file;
            file.dzfname.update(file.inputFld.value = file.name);
        }
        
        var updateDropZone = function(self) {
            var files = self.fileInputEl.dom.files;
            if (files.length == 0) {
                uPanel.down('#dz_' + self.itemId).update('');
                uPanel.uploadButt.disable();
            } else {
                processFName(self, files[0]);
                uPanel.uploadButt.enable();
            }
        };
        var helpTxt = fileUploadDisabled ? 'File upload is unavailable in Internet Explorer versions 9 and earlier' : '';
        if (SeqView.jsonType == 'JSONP') {
            fileUploadDisabled = true;
            helpTxt = 'Local files uploading is currently unavailable for X-domain/IE configuration';
        }
        var diff_parts = {
            rid: [
                {xtype:'displayfield', value: 'Please enter NCBI BLAST request ticket (RID) then press Upload to add new alignment track.',
                     width:'95%',  style: {fontSize: '122%'} },
                {xtype:'textfield', hideLabel:true, emptyText:'Please enter Blast RID', name: 'rid', width:'99%',
                     validator: manageButt},
                {xtype:'checkbox', name: 'find_comp', height:15, boxLabel:'Link related hits together',  checked:true, style: {fontSize: '100%'}},
                {xtype:'displayfield', value: 'BLAST returns separate alignments for each query, and these separate alignments can further be ordered into sets offering consistent non-overlapping query and subject coverage.  The sequence viewer offers the ability to evaluate the original BLAST hits on-the-fly and link together alignments that meet a strict definition of non-overlapping query and subject coverage.',
                     height: 65}
            ],
            text: [
                {xtype:'displayfield', value: 'Please paste you text annotations then press Upload to add new track.',
                     width:'95%',  style:{fontSize: '122%'} },
                {xtype:'textarea', hideLabel:true, emptyText:'Please paste your text here', name: 'data', width:'99%', height:135,
                     validator: manageButt}
            ],
            url: [
                {xtype:'displayfield', value: 'Please specify the WEB URL to download the input file then press Upload to add new track(s).',
                    width:'95%',  style:{fontSize: '122%'} },
                {xtype:'textfield', hideLabel:true, emptyText:'Please enter URL', name: 'dataURL', width:'99%', validator: manageButt}
            ],
            BAMfiles: (helpTxt)
                ? [{xtype:'displayfield', value: helpTxt,
                    width:'95%',  style: {fontSize: '122%', color: 'red'}}]
                : [{xtype:'displayfield', value: 'Please specify or drop input files then press Upload to add new track(s).',
                     width:'95%',  style: {fontSize: '122%'}},
                   {xtype:'panel', layout: 'form',  border: false, bodyStyle: {background: 'inherit'},
                        items: [
                            {xtype: 'filefield',  itemId: 'file', fieldLabel: 'File to upload', listeners: {change: updateDropZone}},
                            {xtype: 'filefield',  itemId: 'auxiliary_file', fieldLabel: 'Index file', listeners: {change: updateDropZone}}
                        ]
                   },
                   {xtype: 'fieldset', title: 'Drag and drop file here', id: 'sv_dropzone', name: 'bam', height: 60, labelWidth: 1, width: '100%',
                        items: [{xtype: 'displayfield', itemId: 'dz_file', value: ''},
                                {xtype: 'displayfield', itemId: 'dz_auxiliary_file', value: ''}]}
                ],
            FASTAfile: (helpTxt)
                ? [{xtype:'displayfield', value: helpTxt,
                    width:'95%',  style: {fontSize: '122%', color: 'red'}}]
                : [{xtype:'displayfield', value: 'Please specify or drop input MUSCLE/FASTA file then press Upload to add new track(s).',
                     width:'95%',  style: {fontSize: '122%'}},
                   {xtype:'hidden', name:'file_format', value: 'align'},
                   {xtype:'panel', layout: 'form',  border: false, bodyStyle: {background: 'inherit'},
                        items: [
                            {xtype: 'filefield',  itemId: 'file', fieldLabel: 'File to upload', listeners: {change: updateDropZone}}
                        ]
                   },
                   {xtype: 'fieldset', title: 'Drag and drop file here', id: 'sv_dropzone', name: 'align', height: 60, labelWidth: 1, width: '100%',
                        items: [{xtype: 'displayfield', itemId: 'dz_file', value: ''}]}
                ],

            file: (helpTxt)
                ? [{xtype:'displayfield', value: helpTxt,
                    width:'95%',  style: {fontSize: '122%', color: 'red'}}]
                : [{xtype:'displayfield', value: 'Please specify or drop an input file then press Upload to add new track(s).',
                     width:'95%',  style: {fontSize: '122%'}},
                   {xtype:'panel', layout: 'form',  border: false, bodyStyle: {background: 'inherit'},
                        items: [
                            {xtype: 'filefield',  itemId: 'file', fieldLabel: 'File to upload', listeners: {change: updateDropZone}},
                            {xtype: 'combo', triggerAction:'all', width:120, name:'file_format', fieldLabel:'File format', mode:'local',
                                store:['auto detect', 'asn binary', 'asn text', 'bed', 'bed 15', 'fasta', 'gff2', 'gff3', 'gtf/gff', 'gvf', 'vcf', 'wiggle', 'xml'],
                                allowBlank: false, editable: false, value: 'auto detect'}
                        ]
                   },
                   {xtype: 'fieldset', title: 'Drag and drop file here', id: 'sv_dropzone', height: 60, labelWidth: 1, width: '100%',
                        items: [{xtype: 'displayfield', itemId: 'dz_file', value: ''}]}
                ]
            };
        var panel = this.down('#loaddatapanel');
        var dp = panel.down('#diff_part');
        if (dp) panel.remove(dp);
        panel.insert(0, new Ext.form.Panel({items: diff_parts[type], itemId: 'diff_part', border: false, bodyStyle: {background: 'inherit'}}));
        panel.setTitle(title);
        panel.updateLayout();

        if (!fileUploadDisabled && type.indexOf('file') >= 0) {   
            var dzone = document.getElementById('sv_dropzone');
            dzone.ondragover = function(){ return false; };
            dzone.ondragleave = function(){ return false; };
            var dtfile = uPanel.down('#file'),
                idxfile = uPanel.down('#auxiliary_file');
            dzone.ondrop = function(e) {
                e.preventDefault();
                var files = e.dataTransfer.files;
                if (type == 'BAMfiles') {
                    Ext.each(files, function(f, idx) {
                        processFName(/.(BAM)$/i.test(f.name) ? dtfile : idxfile, f);
                    }, uPanel, true);
                    if (!uPanel.dzfiles || !uPanel.dzfiles.file || !uPanel.dzfiles.auxiliary_file) return;
                }
                else processFName(dtfile, files[0]);
                uPanel.uploadButt.enable();
                uPanel.cfgPanel.sviewApp.pingClick('9-5-DnD');
                return false;
            };
        }
    }
});
/*  $Id: tm_dialog.js 38325 2017-04-25 22:08:59Z borodine $
*/

function encodeKeyVal(s, key, val, f_escape) {
    if (val && val.length > 0) {
        s += (s.length > 0 ? ',' : '') + key + ':' + (f_escape ?  NCBIGBUtils.escapeTrackName(val) : val);
    }
    return s;
}


SeqView.TM.processTracksInfo = function(cgi_tracks) {
    var ret = [];
    var search_category = function(categories, cat_name) {
        var scat = null;
        Ext.each(categories,function(cat) {
            if (cat.name == cat_name) {
                scat = cat;
                return false;
            }
        });
        return scat;
    }
    var uu_id = 0;
    var next_id = function() {
        return uu_id++;
    }
    Ext.each(cgi_tracks,function(ctrack) {
        var cat = null;
        ctrack.uuid = next_id();
        if( !ctrack.display_name ){
            if( !ctrack.name || ctrack.name.length == 0 ){
                ctrack.name = "track_" + ctrack.uuid;
                ctrack.display_name = "Track " + ctrack.uuid;
            } else {
                ctrack.display_name = "Track " + ctrack.name;
            }
        } else {
            if( !ctrack.name || ctrack.name.length == 0 ){
                ctrack.name = "track_" + ctrack.display_name;
            }
        }
        if (!ctrack.order && ctrack.order !== 0)
            ctrack.order = 100;
        if(ctrack.category) {
            cat = search_category(ret,ctrack.category.name)
            if (!cat) {
                cat = Ext.clone(ctrack.category);
                cat.subcategories = [];
                ret.push(cat);
            }
        } else {
            cat = search_category(ret,"other")
            if (!cat) {
                cat = {name:"other", display_name:"Other", help: "Other tracks"};
                cat.subcategories = [];
                ret.push(cat);
            }
        }
        var subcat = null;
        if (ctrack.subcategory) {
            subcat = search_category(cat.subcategories,ctrack.subcategory.name)
            if (!subcat) {
                subcat = Ext.clone(ctrack.subcategory);
                subcat.tracks = [];
                cat.subcategories.push(subcat);
            }
        } else {
            subcat = search_category(cat.subcategories,"other")
            if (!subcat) {
                subcat = {name:"other", display_name:"Other", help: "Other tracks"};
                subcat.tracks = [];
                cat.subcategories.push(subcat);
            }
        }
        subcat.tracks.push(ctrack);
    });
    ret.sort(function(a, b) { return a.order - b.order; });
    return ret;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//
SeqView.TM.trackToString = function(track, for_group) {

    var str = '';

    str = encodeKeyVal(str, 'key', track.key);
    if (for_group) {
        str = encodeKeyVal(str, 'group_name', track.setting_group);
    } else {
        str = encodeKeyVal(str, 'name', track.name, true);
        str = encodeKeyVal(str, 'display_name', track.display_name, true);
        str = encodeKeyVal(str, 'id', track.id, true);
        str = encodeKeyVal(str, 'data_key', track.data_key, true);
        str = encodeKeyVal(str, 'subkey', track.subkey);
        str = encodeKeyVal(str, 'category', track.cat_name);
        str = encodeKeyVal(str, 'subcategory', track.subcat_name);
        str = encodeKeyVal(str, 'dbname', track.dbname, true);
        str = encodeKeyVal(str, 'uid', track.uId, true);
        Ext.each(track.subTracks, function(subtrack, idx) {
            str += idx == 0 ? ',subtracks:' : '|';
            str += NCBIGBUtils.escapeTrackName(subtrack);
        });
        Ext.each(track.annots, function(annot, idx) {
            str += idx == 0 ? ',annots:' : '|';
            str += NCBIGBUtils.escapeTrackName(annot);
        });
        Ext.each(track.comments, function(comment, idx) {
            str += idx == 0 ? ',comments:' : '|';
            str += NCBIGBUtils.escapeTrackName(comment.label) + '|' + comment.pos_str;
        });
        Ext.each(track.highlights, function(hl, idx) {
            str += idx == 0 ? ',highlights:' : '|';
            str += NCBIGBUtils.escapeTrackName(hl);
        });
        if (track.filter && track.filter.length > 0) {
            // Do not escape the bar '|'
            str += ',filter:' + NCBIGBUtils.escapeTrackName(track.filter, /([\]\[\\,:=&;"#%])/);
        }
    }
    str += SeqView.TM.getTrackDisplayOptions(track);
    if (track.highlights_color) str += ',highlights_color:' + track.highlights_color;
    if (track.is_private) str += ',is_private:true';
    if (track.stored_scale) str += ',stored_scale:' + track.stored_scale;
    return str;
}

SeqView.TM.getTrackDisplayOptions = function(track) {
    var str = '';
    Ext.each(['choice_list', 'check_boxes', 'text_boxes', 'hidden_settings'], function() {
        Ext.each(track[this], function() {
            str += ',' + this.name + ':' + (this.hasOwnProperty('curr_value') ? this.curr_value : this.value);
        });
    });
    return str;
}

SeqView.TM.tracksArrayToString = function(tracks, shown, order, extraOpt) {
    var str = '';
    if (typeof extraOpt != 'string') extraOpt = '';
    Ext.each(tracks, function(track) {
        str += '[' + SeqView.TM.trackToString(track, null);
        if (shown) str += ',shown:' + (track.shown ? 'true' : 'false');
        if (order) str += ',order:' + track.order;
        str += extraOpt + ']';
    });
    return str;
}

// returns an array of strings for parallel rendering
SeqView.TM.tracksToArrayOfStrings = function(tracks, extraOpt){
    if (!tracks.length) return [];
    var newGroup = curGroup = -1;
    var trx = [''];

    for (var i = 0; i < tracks.length; i++) {
        if (!tracks[i].shown) continue;
        if (curGroup != tracks[i].render_group) {
            curGroup = tracks[i].render_group;
            if (++newGroup) trx[newGroup] = '[key:no ruler]';
        }
        var trk = tracks[i];
        trk.rg = newGroup;
        trx[newGroup] += '[' + SeqView.TM.trackToString(trk, null) + extraOpt + ']';
        if (trk.key == 'graph_overlay' && trk.subTracks) {
            trk.subTracks.forEach(function(uId, idx){
                if (uId.substr(-7) == '_hidden') return;
                for (var j = 0; j < tracks.length; j++) {
                    if (tracks[j].id != trk.legend[idx].id) continue;
                    if (trk.render_group != tracks[j].render_group || !tracks[j].shown)
                        trx[newGroup] += '[' + SeqView.TM.trackToString(tracks[j], null)+ ', shown:false' + extraOpt + ']';
                    break;
                }
            });
        }
    }
    if (SeqView.TM.renderStat) {
        tracks.forEach(function(t) {
            console.log(t.display_name + ', render_group: ' + t.render_group + ', rg# ' + t.rg + ', order: ' + t.order);
        });
    }
    return trx;
}


SeqView.TM.generateTracksString = function(categories){
    categories = categories || SeqView.TM.categories;
    if (!categories) return null;

    var act_tracks = SeqView.TM.collectActiveTracks(categories);
    return SeqView.TM.tracksArrayToString(act_tracks, true);
}


SeqView.TM.generateTracksStringForConfig = function(categories, filter_shown) {
    categories = categories || SeqView.TM.categories;
    if (!categories) return null;

    var str = '';
    Ext.each(categories, function(cat) {
        Ext.each(cat.subcategories, function(subcat) {
            Ext.each(subcat.tracks, function(track) {
                if (track.shown || !filter_shown) {
                    str += '[';
                    str += SeqView.TM.trackToString(track);
                    str += ',order:' + track.order;
                    str += ',shown:' + (track.shown ? 'true' : 'false');
                    str += ']';
                }
            });
        });
    });
    return str;
}


SeqView.TM.applyGroupSettings = function(categories, group_track) {
    Ext.each(categories, function(cat) {
        Ext.each(cat.subcategories, function(subcat) {
            Ext.each(subcat.tracks, function(track) {
                if (track.key === group_track.key &&
                    track.setting_group === group_track.setting_group) {
                    track.check_boxes = group_track.check_boxes;
                    track.choice_list = group_track.choice_list;
                    if (group_track.text_boxes) track.text_boxes = group_track.text_boxes;
                }
            });
        });
    });
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
SeqView.TM.getGroupSettingsString = function(group_track) {
    var str = '';
    if (group_track) {
        str += '[';
        str += SeqView.TM.trackToString(group_track, true);
        str += ']';
    }
    return str;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
SeqView.TM.Common = {
    updateSeqViewApp: function(categories, app, group_track) {
        var group_settings;
        app = app || this.sviewApp;
        var user_data = app.m_Key;
        if (this.getEl) this.getEl().mask('Configuring...');
        var seqConfig = this.seqConfig || app.m_Config;
        categories = categories || SeqView.TM.processTracksInfo(seqConfig.TrackConfig);
        group_track = group_track || this.group_track || null;
        if (group_track) {
            SeqView.TM.applyGroupSettings(categories, group_track);
            group_settings = SeqView.TM.getGroupSettingsString(group_track);
        }
        seqConfig.save({
            tracks: (typeof categories == 'string') ? categories : SeqView.TM.generateTracksStringForConfig(categories),
            group_settings: group_settings,
            callback: {
                success: function(res) {
                    if (!user_data || user_data.length == 0) user_data = null;
                    var tracks = SeqView.TM.generateTracksString(categories);
                    if (!tracks || tracks.length == 0) tracks = null;
                    if (this.defaultTracks) app.forEachView(function(v){ delete v.m_showAllGenes; });
                    if (this.cleanCfg) this.cleanCfg();
                    app.updateConfig(seqConfig, user_data, tracks);
                },
                failure: function(res) {
                    var fromcgi = Ext.decode((typeof res === 'object') ? res.responseText : res);
                    var msg = fromcgi.success === false ? fromcgi.msg : res.responseText;

                    Ext.MessageBox.show({title: 'Configuration saving error',msg: msg,
                                 buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                    if (this.cleanCfg) this.cleanCfg();
                },
                scope: this
            }
        });
    },

    createMainPanel: function(categories) {
        var removeButton = function() {
            var agp = this.tracksPanel.getActiveGridPanel();
            var win = new Ext.Window({title: 'Confirmation', modal: true, labelWidth: 1,
                app: this.sviewApp,
                layout:'form', width:460, buttonAlign: 'center',
                buttons:[
                    {text: 'Continue', handler: function() {
                        win.close();
                        this.removeUUDTracks(); }, scope: this},
                    {text: 'Cancel', handler: function() { win.close(); }}]
            });
            win.add({xtype: 'displayfield', value: 'The following track(s) will be permanently removed from NCBI Data Storage after "Configure" process is submitted:',
                style: { 'font-size': '13px' }});

            Ext.each(agp.getSelection(), function(t, idx){
                var trk = t.getData();
                if (trk.id.charAt(0) != 'U') 
                    agp.getSelectionModel().deselect(t);
                else
                    win.add({xtype: 'displayfield', value: (idx + 1) + '. ' + trk.display_name,
                         style: { 'font-weight': 'bold', 'font-size': '13px',
                        'background-color': (idx%2 ? '' : 'gainsboro')}});
            });
            win.show();
        }

        this.categories = categories;
        this.tracksPanel = new SeqView.TM.MainPanel({
            sviewApp: this.sviewApp,
            categories: this.categories
        });

        this.uploadPanel = new SeqView.TM.UploadPanel({
            cfgPanel: this,
            dirty: false

        });

        var activeTab = this.activeTab || 0;
        var panel = new Ext.TabPanel({
            activeTab: activeTab,
            buttonAlign: 'left',
            items: [this.tracksPanel, this.uploadPanel],
            listeners: {
                beforetabchange: function(p, newtab, curtab) {
                    var isTP = newtab === this.tracksPanel;
                    Ext.getCmp('sv-delete_button_id')[isTP ? 'show' : 'hide']();
                    if (isTP && newtab.dirty) {
                        this._loadTracks(false);
                    }
                },
                scope: this
            },
            buttons: [
                {text: 'Remove track(s)', scope: this, icon: SeqView.base_url + 'images/delete.png',
                tooltip: 'Remove track(s) from NCBI data storage',
                disabled: true, id: 'sv-delete_button_id',
                handler: removeButton},
                '->',
                {text: 'Configure Group', scope: this,
                hidden: true, id: 'sv-config_group_button_id',
                handler: function(b) {
                    if (this.cfgWindow) {
                        if (this.sviewApp) this.group_track = b.track;
                        this.reConfigure(false);
                    }
                }},
                {text: 'Configure', id: 'sv-configure_button_id', scope: this, handler: function(b) {
                    this.sviewApp.cfgEvent = true;
                    if (!this.uploadPanel.uploadButt.disabled) {
                        var msgBox = Ext.MessageBox.show({title:'External Data Upload',
                            msg: 'Upload external data?',
                            buttons:Ext.MessageBox.YESNO,
                            icon:Ext.MessageBox.QUESTION,
                            fn: function(buttonId) {
                                this.reConfigure(buttonId == 'yes');
                                delete msgBox;
                            },
                            scope: this
                        });
                    }
                    else this.reConfigure(false);
                }},
                {text: 'Load Defaults', id: 'sv-load_default_button_id', scope: this, handler: function(b) {
                    this._loadTracks(this.defaultTracks = true);
                }},
                {text: 'Cancel', scope: this, handler: function() {
                    if (this.cleanCfg()) return;
                    this.seqConfig.TrackConfig = Ext.clone(this.sviewApp.m_Config.TrackConfig)
                    this.seqConfig.Options = Ext.clone(this.sviewApp.m_Config.Options);
                    this.categories = SeqView.TM.processTracksInfo(this.seqConfig.TrackConfig);
                    this.tracksPanel.reload(this.categories);
                }}
            ]
        });
        return panel;
    },

    removeUUDTracks: function() {
        var agp = this.tracksPanel.getActiveGridPanel(),
            aStore = this.tracksPanel.activeTracksTab.gridPanel.store,
            sStore = this.tracksPanel.searchTracksTab.gridPanel.store;
        var sm = agp.getSelectionModel();
        var selection = sm.getSelection();
        this.removeList = selection.concat(this.removeList || []);
        Ext.each(selection, function() {
            var uuid = this.getData().uuid;
            aStore.removeAt(aStore.find('uuid', uuid));
            sStore.removeAt(sStore.find('uuid', uuid));
            if (sm.lastSelected == this && !sm.selectNext()) sm.selectPrevious();
            agp.store.remove(this);
        });
    },
    // private
    _loadTracks: function(default_cfg) {
        this.getEl().mask('Loading...');
        if (default_cfg)
            this.sviewApp.m_Key = this.sviewApp.defaultConfig.m_Key;
        this.seqConfig.load({
            defaultcfg: default_cfg,
            user_data: this.sviewApp.m_Key,
            callback: {
                success: function(config) {
                    this.categories = SeqView.TM.processTracksInfo(config.TrackConfig);
                    this.tracksPanel.reload(this.categories);
                    this.getEl().unmask();
                },
                failure: function(res) {
                    var fromcgi = Ext.decode((typeof res === 'object') ? res.responseText : res);
                    var msg = fromcgi.success === false ? fromcgi.msg : res.responseText;

                    Ext.MessageBox.show({title: 'Configuration loading error',msg: msg,
                                         buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                    this.getEl().unmask();
                },
                scope: this
            },
            forcereload: true
        });
    }
}

SeqView.TM.ShowConfigDialog = function(sview_app, activeTab) {
    sview_app.resizeIFrame(600);
    var seqconfig = new SeqView.Config(sview_app);
    seqconfig.TrackConfig = Ext.clone(sview_app.m_Config.TrackConfig);
    seqconfig.Options = Ext.clone(sview_app.m_Config.Options);
    var tm = new Ext.Window({
        title: 'Configure Page',
        app: sview_app,
//        modal: true,
        layout:'fit',
        width:750,
        height:550,
        plain: true,
        cls: 'SeqViewerApp'
    });
    var cfgPanel = Ext.create('SeqView.TM.ConfigPanel', {
       cfgWindow: tm,
       accession: sview_app.m_Config.SeqInfo.id,
       seqConfig: seqconfig,
       sviewApp: sview_app,
       activeTab: activeTab || 0
    });
    tm.add(cfgPanel);
    tm.on('close', function(p) {
        Ext.getBody().unmask();
        sview_app.resizeIFrame();
        sview_app.m_DialogShown = false;
        sview_app.fireEvent('configuration_panel', sview_app, sview_app.cfgEvent ?  'close' : 'cancel');
        delete sview_app.cfgEvent;
    });
    sview_app.fireEvent('configuration_panel', sview_app, 'open');
    sview_app.m_DialogShown = true;

    if (sview_app.m_iFrame) Ext.defer(tm.show, 500, tm); else tm.show();
    Ext.getBody().mask();
};


SeqView.TM.modifyTrackDetails = function(gview, track) {
    var app = gview.m_App;
    app.resizeIFrame(600);
    var tsWindow = new Ext.Window({
        title: track.display_name,
        app: app,
        modal: true,
        width:550,
        cls: 'SeqViewerApp'
    });

    var detailPanel = new SeqView.TM.TrackDetail({
        region: 'north',
        height: 250,
        minSize: 100/*,
        changedGenesMode: false,
        gview: gview*/
    });

    var applySettings = function(button) {
        if (track.category.name == 'Genes') gview.m_App.m_Config.cleanGeneMode(track);

        SeqView.TM.Common.updateSeqViewApp(null, app, button.group_tracks);
        tsWindow.close();
    }

    tsWindow.add(detailPanel);
    var gname = track.setting_group;
    if (gname) {
        var display_limit = 10;
        if (gname.length > display_limit) gname = group.substr(0, display_limit) + "...";
    }
    var panel = tsWindow.add({buttons: [
        {text: "Apply to Group \'" + gname + "'", scope: this, handler: applySettings,
         hidden: !(gname), group_track: track},
        {text: 'Apply',  scope: this, handler: applySettings},
        {text: 'Cancel', scope: this, handler: function() {tsWindow.close();}}]});

    SeqView.TM.updateTrackDetails(detailPanel, track);
    tsWindow.on('render', function(){
        panel.getDockedItems()[0].insert(0, [{ xtype: 'component',
            html: '<a onclick="SeqView.showHelpDlg(\'legends/#' + track.legend_text + '\');" href="#">Track legend</a>'}, '->']);
    });
    tsWindow.on('close', function() {
        app.resizeIFrame();
    });
    if (app.m_iFrame) Ext.defer(tsWindow.show, 500, tsWindow); else tsWindow.show();
};


Ext.define('SeqView.TM.ConfigPanel', {
    extend: 'Ext.Panel',
    initComponent : function(){
        var categories = SeqView.TM.processTracksInfo(this.seqConfig.TrackConfig);
        Ext.apply(this, {
            border: false,
            height:  515,
            layout: 'fit',
            items: this.createMainPanel(categories)
        });
        this.callParent(arguments);
    },
    cleanCfg: function() {
        if (this.getEl().isMasked()) this.getEl().unmask();
        delete this.defaultTracks;
        if (!this.cfgWindow) return false;
        this.cfgWindow.close();
        return true;
    },
    reConfigure: function(checkDirty) {
        if (this.sviewApp) {
            if (checkDirty && !this.uploadPanel.uploadButt.disabled) {
                this.uploadPanel.submitData(this.down('#loaddatapanel').getForm(),
                                            this.cfgProcess, this);
            } else {
                this.cfgProcess();
            }
        } else {
            this.cleanCfg();
        }
    },
    cfgProcess: function() {
//        this.sviewApp.forEachView(function(v) { v.setGeneMode(false); }); 
        if (this.tracksPanel.dirty) {
            this.getEl().mask('Loading...');
            this.seqConfig.load({
                user_data: this.sviewApp.m_Key,
                callback: {
                    success: function(config) {
                        var categories = SeqView.TM.processTracksInfo(config.TrackConfig);
                        this.updateSeqViewApp(categories);
                    },
                    failure: function(res) {
                    var fromcgi = Ext.decode((typeof res === 'object') ? res.responseText : res);
                        var msg = fromcgi.success === false ? fromcgi.msg : res.responseText;
                        delete this.sviewApp.cfgEvent;
                        Ext.MessageBox.show({title: 'Configuration loading error',msg: msg,
                                           buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                    },
                    scope: this
                },
                forcereload: true
            });
        } else {
            var seqConfig = this.sviewApp.m_Config;
            var trkCfg = this.seqConfig.TrackConfig;
            Ext.each(this.removeList, function(){
                var uuid = this.data.uuid;
                for (var j = trkCfg.length - 1; j >= 0; j--) {
                    if (uuid == trkCfg[j].uuid) { trkCfg.splice(j, 1); break; }
                }
                seqConfig.removeUserTrack(this.data.uuid);
            });
            var agp = this.tracksPanel.getActiveGridPanel();
            if (agp.syncModifiedTracks() && this.tracksPanel.searchTracksTab.gridPanel == agp) {
                // if the active is Search Tracks tab and there were changes it overwrites changes made in other tabs
                var tcfg = this.tracksPanel.sviewApp.m_Config.TrackConfig;
                this.categories = SeqView.TM.processTracksInfo(tcfg);
                this.seqConfig.TrackConfig = Ext.clone(tcfg);
            }
            this.updateSeqViewApp(this.categories);
            if (this.sviewApp.m_PermConfId) this.tracksPanel.reload(this.categories);
        }
    }
});

Ext.apply(SeqView.TM.ConfigPanel.prototype, SeqView.TM.Common);

SeqView.TM.ShowConfigPanel = function(sview_app, activeTab) {
    var seqconfig = new SeqView.Config(sview_app);
    seqconfig.TrackConfig = Ext.clone(sview_app.m_Config.TrackConfig);
    seqconfig.Options = Ext.clone(sview_app.m_Config.Options);
    var confpanel = Ext.create('SeqView.TM.ConfigPanel', {
        accession: sview_app.m_Config.SeqInfo.id,
        seqConfig: seqconfig,
        sviewApp: sview_app
    });
    var tm = new Ext.Panel({
        title: 'Configure Page',
        layout: 'fit',
        border: false,
        renderTo: sview_app.m_PermConfId,
        items: confpanel
    });

    var el = Ext.get( sview_app.m_PermConfId );
    if( el ){
        el.set({ 'confpanel_id': confpanel.id });
    }
    tm.show();
};

SeqView.TM.modifyLegendSettings = function(gview, track, idx) {
    var app = gview.m_App;
    var legends = track.legends;
    if (!legends[0].name)
        Ext.each(legends, function() {
            var split = this.label.split(': ');
            if (split.length > 2) split[1] = this.label.substr(split[0].length + 2);
            this.gStyle = split[0];
            this.name = split[1];
            this.hidden = !this.bounds;
        });
    app.resizeIFrame(300);
    var lgWindow = new Ext.Window({
        title: track.display_name,
        app: app,
        modal: true,
        width:550,
        listeners: {
            close: function() { app.resizeIFrame(); },
            show: function(){ updateSettings(); }
        }
    });
    var applySettings = function(button) {
        var doReload = false;
        var trackConfig = app.m_Config.TrackConfig;
        var subTracks = [];
        Ext.each(track.legends, function(lg) {
            if (!lg.bounds !== lg.hidden) {
                doReload = true;
                delete lg.bounds;
            }
            subTracks.push(trackConfig[lg.idx].uId + (lg.hidden ? '_hidden' : ''));
            Ext.each(trackConfig[lg.idx].hidden_settings, function() {
                if (lg[this.name] != this.value) {
                    doReload = true;
                    this.value = lg[this.name];
                }
            });
            if (lg.label.split(': ')[0] != lg.gStyle) {
                doReload = true;
                Ext.each(trackConfig[lg.idx].choice_list, function() {
                    if (this.name != 'style') return true;
                    this.curr_value = lg.gStyle;
                    return false;
                });
            }
        });
        if (doReload) {
            trackConfig[track.idx].subTracks = subTracks;
            SeqView.TM.Common.updateSeqViewApp(null, app);
        }
        lgWindow.close();
    }

    var trackStore = [];
    Ext.each(legends, function() { trackStore.push(this.name); });
    var trackPicker = new Ext.form.ComboBox({ triggerAction: 'all', width: 530, name: 'trackPicker', mode:'local',
                       store: trackStore, allowBlank: false, editable: false, value: legends[idx].name,
                       listeners: {select: updateSettings}});

    lgWindow.add(trackPicker);

    var settingsPanel =  new Ext.Panel({
        bodyStyle: {background: 'inherit'},
        padding: '0 3 0 3'
    });
    lgWindow.add(settingsPanel);

    var colorPicker = {xtype: 'container', layout: 'hbox',  hei_ght: 22, padding: '3', // bodyStyle: {background: 'inherit'},
        items: [
           {xtype:'displayfield',  labelAlign: 'right', fieldLabel: 'Color' },
           {xtype:'component', itemId: 'colorBox', width:22, height: 22},
           {xtype:'button',  height: 22,
                menu: { xtype: 'colormenu',
                    handler: function(picker, choice) {
                        legends[idx].color = '#' + choice;
                        settingsPanel.down('#colorBox').el.applyStyles('background-color: ' + legends[idx].color + ';');
                    }
                }
            }]
    };
    var slidePicker = {xtype: 'container', layout: 'hbox',  height: 22, //bodyStyle: {background: 'inherit'},
        items: [
            {xtype: 'slider', itemId: 'opacity', labelAlign: 'right', fieldLabel: 'Opacity', width: '50%', increment: 5,
                tipText: function(thumb) {
                    legends[idx].opacity = thumb.value;
                    settingsPanel.down('#colorBox').el.applyStyles('opacity:' + (thumb.value/100) + ';');
                    return String(thumb.value) + '%';
            }}
       ]
    };
    var stylePicker = new Ext.form.RadioGroup({fieldLabel: 'Graph style', columns: 1, labelAlign: 'right',
            items: [
                {boxLabel: 'Histogram', name: 'gStyle', inputValue: 'histogram'},
                {boxLabel: 'Line graph', name: 'gStyle', inputValue: 'line graph'}
            ],
            listeners: {change: function(sp, val) { legends[idx].gStyle = val.gStyle; }}
     });

    var trackHider = new Ext.form.Checkbox({ boxLabel: 'Hide track', padding: '5',
        listeners: {change: function(o, val) { settingsPanel[(legends[idx].hidden = val) ? 'disable' : 'enable'](); }}
    });
    lgWindow.add(trackHider);

    function updateSettings(tPicker) {
        if (tPicker) {
            var rv = tPicker.getRawValue();
            for (var i = 0; i < legends.length; i++) {
                if (legends[i].name != rv) continue;
                idx = i;
                break;
            }
        }
        var legend = legends[idx];
        settingsPanel.down('#colorBox').el.applyStyles('background-color:' + legend.color + ';opacity:' + (legend.opacity/100) + ';');
        settingsPanel.down('#opacity').setValue(parseInt(legend.opacity));
        stylePicker.setValue({gStyle: legend.gStyle});
        trackHider.setValue(legend.hidden);
    }

    settingsPanel.add(colorPicker, slidePicker, stylePicker);

    lgWindow.add({buttons: [
        {text: 'Apply',  scope: this, handler: applySettings},
        {text: 'Cancel', scope: this, handler: function() {lgWindow.close();}}]});

    if (app.m_iFrame) Ext.defer(lgWindow.show, 500, lgWindow); else lgWindow.show();
};


/*  $Id: tm_grids.js 36981 2016-11-21 21:02:47Z borodine $
*/

Ext.define('SeqView.TrackList', {
    extend: 'Ext.data.Model',
    proxy: {type: 'memory'},
    fields: [
        'display_name',
        'help',
        { name: 'order', type: 'int' },
        { name: 'uuid', type: 'int' },
        { name: 'shown', type: 'bool' },
        { name: 'subcat', mapping: function(data)
            {return data.subcat || data.category.display_name;} }
    ]
});

Ext.define('SeqView.TM.GridBase', {extend: 'Ext.grid.GridPanel',
    border: false,
    itemId: 'gridbase',
    multiSelect: true,
    overflowY: 'auto',
    region: 'center',
/*    listeners: {
        itemmouseenter: function(view, record, item) {
            Ext.fly(item).set({'data-qtip': record.get('help')});
        }
    },*/
    initComponent: function() {
        if (this.subcategories && this.subcategories.length > 1) {
        // groupping tracks panel
            var tracks = [];
            Ext.each(this.subcategories, function(subcat) {
                Ext.each(subcat.tracks, function(track) {
                    track.subcat = subcat.display_name;
                    tracks.push(track);
                });
            });
            this.tracks = tracks;
            this.groupGrid = true;
        }
//        this.viewConfig = { listeners: { render: this.createTooltip } };
        this.store = this.store || Ext.create('Ext.data.Store', {
            autoDestroy: true,
            model: 'SeqView.TrackList',
            data: this.tracks
        });
        if (this.groupGrid) {
            this.store.setGroupField ('subcat');
            this.features = [{
                ftype: 'grouping',
                groupHeaderTpl: 'Category: {name} ({[values.rows.length]} Item{[values.rows.length > 1 ? "s" : ""]})',
                hideGroupedHeader: true,
                enableGroupingMenu: false
            }];
        }
        
        var sortable = this.hasOwnProperty('sortable') ? this.sortable : false;

        if (!this.columns) {
            this.columns = [
                {header: 'Active', width: 45, xtype: 'checkcolumn', sortable: sortable, dataIndex: 'shown', menuDisabled: true,
                    listeners: {checkchange: this.checkChange},
                    editor: { xtype: 'checkbox', cls: 'x-grid-checkheader-editor'}},
                {header: "Track name", sortable: sortable, dataIndex: 'display_name', flex: true, editable: false, menuDisabled: true}];
                if (this.searchTab) {
                    this.columns.unshift({header: "Category", width: 70, sortable: true, dataIndex: 'subcat', editable: false, menuDisabled: true});
                    this.columns.push({header: "Description", dataIndex: 'help', flex: true, editable: false, menuDisabled: true});
                }
        }
        if (this.groupGrid) this.columns.push({header: 'Category', dataIndex: 'subcat', sortable: sortable, width: 70});
        this.callParent(arguments);
    },
    
    checkChange: function(obj, idx, val) {
        var gp = this.up('#gridbase'),
            sel = gp.getSelection(),
            uuid = gp.getStore().getAt(idx).getData().uuid;
        if (Ext.each(sel, function(){ return this.getData().uuid != uuid; }) !== true)
            Ext.each(sel, function(){ this.set('shown', val); });
    },
    
    createTooltip: function(view) {
    	view.tip = Ext.create('NCBIGBObject.ToolTip', {
            target: view.el,
            delegate: view.itemSelector,
            pinnable: false,
            header: false,
            adjustWidth: true,
            listeners: {
                beforeshow: function (tip) {
                    var rec = view.getRecord(tip.triggerElement);
                    var tooltip = rec.get('help') || '';
/*                    var ltxt = rec.get('legend_text') || '';
                    if (ltxt) {
                        tooltip = '<a href=\"' + SeqView.getHelpURL() + 'legends/#' + ltxt + 
                                  '\" target=\"_blank\" style=\"color:blue\">Track legend</a>'
                                + '<div>' + tooltip + '</div>';
                    }*/
                    if(tooltip){
                        tip.update(tooltip);
                    } else {
                         tip.on('show', function(){
                             Ext.defer(tip.hide, 10, tip);
                         }, tip, {single: true});
                    }
                }
            }/*,
            onShow: function() {
//                this.callParent(arguments);
                var rec = view.getRecord(this.triggerElement);
                var tooltip = rec.get('help') || '';
                var ltxt = rec.get('legend_text') || '';
                if (ltxt) {
                    tooltip = '<a href=\"' + SeqView.getHelpURL() + 'legends/#' + ltxt + 
                              '\" target=\"_blank\" style=\"color:blue\">Track legend</a>'
                            + '<div>' + tooltip + '</div>';
                }
                if(tooltip) this.update(tooltip);
            }*/
        });
    },

    syncModifiedTracks: function(tracks_config) {
        var mr = this.getStore().getModifiedRecords();
        tracks_config = tracks_config || this.tracks;
        for (var i = 0; i < mr.length; i++) {
            var data = mr[i].data;
            Ext.each(tracks_config, function(tr) {
                if (tr.uuid == data.uuid) { tr.shown = data.shown; return false; }
            });
        }
        return mr.length;
    }
});


SeqView.TM.createSearchTracksTab = function(title, tracks_config, str) {
    var trackFilter = function(rec) {
        if (!this.baseParams.regexp) return true;
        try {
            if (rec.data.help.search(this.baseParams.regexp) >= 0) return true;
            if (rec.data.display_name.search(this.baseParams.regexp) >= 0) return true;
        } catch(err) { console.log(err.message); }        
        return false;
    };
    var pageSize = 10; 
    var panel = new SeqView.TM.GridBase({
        title: (title) ? title : null,
        border: false,
        searchTab: true,
        sortable: true,
        stripeRows: true,
        tracks: tracks_config,
        viewConfig: {forceFit:true, deferEmptyText:false, emptyText:'<div align="center">No Search Results To Display</div>'}
    });
    var store = panel.store;
    store.baseParams = {query: str || '', limit: pageSize, filterBy: trackFilter, regexp: str ? globStringToRegex(str) : ''};
    if (str) store.filterBy(trackFilter, store);
    if (title) panel.iconCls = (store.getCount()) ? 'xsv-search-results' : '';
    return panel;
};


/*  $Id: sviewapp.js 38378 2017-05-01 18:05:34Z borodine $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 *
 * File Description:
 *
 */

SeqView.fireEvent = function(eName, elemID){
    var elem = document.getElementById(elemID);
    if (document.createEvent) {
        var e = document.createEvent('HTMLEvents');
        e.initEvent(eName, false, false);
        elem.dispatchEvent(e);
    } else elem.fireEvent(eName);
};


/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Config
/********************************************************************/

SeqView.Config = (function() {

    function constructor(app) {
        this.m_App = app;
    }

    return constructor;

}) ();

SeqView.Config.prototype = {
    save: function(cfg) {
        var succ_fn = cfg && cfg.callback && cfg.callback.success ? cfg.callback.success : function(seqconfig) {};
        var fail_fn = cfg && cfg.callback && cfg.callback.failure ? cfg.callback.failure : function() {};
        var scope = cfg && cfg.callback && cfg.callback.scope ? cfg.callback.scope : this;

        var url = cfg.config_url || this.m_App.m_CGIs.Config;
        var cookie_name = SeqView.Cookies.UserTracksCookieName;
        if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserTracksCookieNameBase) {
            cookie_name += '-' + this.m_App.m_AppName;
        }
        var user_config_key = cfg.user_config_key || SeqView.UserTracks.get( cookie_name, null );
        if (!cfg.tracks) {
            succ_fn.call(scope, this); // nothing to save - just call success function and exit
            return;
        }
        var params = {saveconfig: true, tracks: cfg.tracks, userconfigkey: user_config_key};
        if (cfg.group_settings) {
            params.group_settings = cfg.group_settings;
        }
        Ext.apply(params, this.visualOptionsUrl());

        var req_num = this.m_App.m_ReqNum || 0;
        this.m_App.m_ReqNum = ++req_num;

        this.m_App.m_timeStamps = {prior2seqconfig: new Date().getTime()};

        this.m_App.AjaxRequest({url:url, context: this, data: params,
            success: function(data, text, res) {
                if (this.m_App.m_ReqNum > req_num) return;
                var from_cgi = SeqView.decode(data);
                if (from_cgi.success === false) {
                    fail_fn.apply(scope, arguments);
                } else {
                    if (!cfg.user_config_key)
                        SeqView.UserTracks.set(cookie_name, from_cgi.userconfigkey);
                    else
                        this.userconfigkey = from_cgi.userconfigkey;
                    succ_fn.call(scope, this);
                }
            },
            error: function(data, text, res) {
                if (this.m_App.m_ReqNum > req_num) return;
                fail_fn.apply(scope, arguments);
            }
        });

    },

    load: function(cfg) {
        cfg = cfg || {};
        var succ_fn = cfg.callback && cfg.callback.success ? cfg.callback.success : function(seqconfig) {};
        var fail_fn = cfg.callback && cfg.callback.failure ? cfg.callback.failure : function() {};
        var scope = cfg.callback && cfg.callback.scope ? cfg.callback.scope : this;

        var url = cfg.config_url || this.m_App.m_CGIs.Config;
        var seqid = cfg.seqid || this.m_App.GI;

        var params = { id: seqid };
        if (this.m_App.defaultConfig.nocache) params.nocache = this.m_App.defaultConfig.nocache;
        // Save initial parameters for Load Defaults, or restore them from saved
        // All other relevant parameters seem to be untouched
        if (cfg.defaultcfg) {
            // Restore
            this.m_App.m_TracksFromURL = this.m_App.defaultConfig.m_TracksFromURL;
            this.m_App.m_ViewContent   = this.m_App.defaultConfig.m_ViewContent;
            if (this.m_App.defaultConfig.uploadedIDs)
                this.m_App.uploadedIDs = this.m_App.defaultConfig.uploadedIDs.slice(0);
        } else {
            // Save
            /* never happens
            if (!this.m_App.defaultConfig) {
                this.m_App.defaultConfig = {
                    'm_TracksFromURL' : this.m_App.m_TracksFromURL,
                    'm_ViewContent' : this.m_App.m_ViewContent
                };
                if (this.m_App.uploadedIDs)
                    this.m_App.defaultConfig.uploadedIDs = this.m_App.uploadedIDs.slice(0);
                if (this.m_App.m_Key)
                    this.m_App.defaultConfig.m_Key = this.m_App.m_Key;
            }*/
        }

        Ext.each(
            [
                ['m_TracksFromURL', 'tracks'],
                ['m_ViewTheme', 'theme'],
                ['m_ViewContent','content'],
                ['m_ItemID', 'itemID'],
                ['m_AppContext', 'app_context'],
                ['m_AssmContext', 'assm_context'],
                ['m_noGuessAssm', 'noguess_assm'],
                ['m_SRZ', 'srz'],
                ['m_BamPath', 'bam_path'],
                ['m_DepthLimit', 'depthlimit']
            ],
            function( pair ) {
                if( this.m_App[pair[0]] ){
                    params[pair[1]] = this.m_App[pair[0]];
                }
            },
            this
        );

        if (cfg.resettracks || (cfg.defaultcfg && !this.m_App.m_TracksFromURL)) {
            delete params.tracks;
            var cookie_name = SeqView.Cookies.UserTracksCookieName;
            if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserTracksCookieNameBase) {
                cookie_name += '-' + this.m_App.m_AppName;
            }
            params.clear = cookie_name;
        }
        else if (!params.tracks) {
            params.tracks = !this.TrackConfig ? '[amend]'
                            : SeqView.TM.generateTracksStringForConfig(SeqView.TM.processTracksInfo(this.TrackConfig), true);
        }
        if (this.m_App.uploadedIDs) {
            Ext.each(this.m_App.uploadedIDs, function(track) { 
                var id = track.id || track.GetTMSId();
                // test for id:id enclosed by combination of commas or square brackets
                // /(\[|,){1}id:123(\]|,){1}/
                var id_regex = new RegExp("(\\[|,){1}id:" + id + "(\\]|,){1}");
                if (id_regex.test(params.tracks)) // track with this id already exists
                    return;
                // any id parameters
                // (\[|,){1}id:\w+(\]|,){1}
                id_regex = /(\[|,){1}id:\w+(\]|,){1}/;
                
                var not_found = true;
                var name = track.GetName();
                var key = track.GetAttr('track_type');
                if (name !== "" && key !== "") {
                    name = 'annots:' + name;
                    key = 'key:' + key;
                    var tr = params.tracks.split(/\[(.*?)\]/).filter(Boolean);
                    for (var i = 0; i < tr.length; i++) {
                        if (id_regex.test(tr[i])) // track has TMS id already
                            continue;
                        var p = tr[i].split(',');
                        if (p.indexOf(name) != -1 && p.indexOf(key) != -1) {
                            tr[i] += ',id:' + id;
                            not_found = false;
                        }
                    }                
                }
                if (not_found)
                    params.tracks += '[id:' + id + ']';
                else {
                    params.tracks = '[' + tr.join('][') + ']';
                }
            });
            delete this.m_App.uploadedIDs;
        }
        // We only use content param from URL ones for the first load.
        // After that we will get it for m_Config.Options.curr_content
        delete this.m_App.m_ViewContent;
        // Same. Get track configuration from TrackConfig
        delete this.m_App.m_TracksFromURL;

        var viewer_context = cfg.viewer_context || this.m_App.m_ViewerContext;
        if( viewer_context ) params.viewer_context = viewer_context;

        if(this.TrackConfig && !cfg.forcereload ) params.notrackconfig = 1;

        if( this.m_App.m_AppName && this.m_App.m_AppName.length > 0 ){
            params.appname = this.m_App.m_AppName;
        }

        if( cfg.defaultcfg !== true ){
            Ext.apply( params, this.visualOptionsUrl() );

            var cookie_name = SeqView.Cookies.UserTracksCookieName;
            if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserTracksCookieNameBase) {
                cookie_name += '-' + this.m_App.m_AppName;
            }

            var user_config_key = cfg.user_config_key || SeqView.UserTracks.get(cookie_name, null);
            if( user_config_key ) params.userconfigkey = user_config_key;
        } else {
            // we need to pass app specific data in case we have one (ISCA, A2A browsers);
            var app_data = SeqView.SessionData.get( SeqView.Cookies.AppDataCookieName, null );
            if( app_data ) params.key = app_data;

            if( this.Options && this.Options.curr_content ){
                params.content = this.Options.curr_content;
            }

            params.nocache = 1;
        }

        var user_data = cfg.user_data || this.m_App.m_Key;
        if (user_data) params.key = user_data;

        var req_num = this.m_App.m_ReqNum || 0;
        this.m_App.m_ReqNum = ++req_num;

        var parKey = params.key;
        var delay = 0;

        if (this.m_App.m_Panel.tmpPanel) this.m_App.m_Panel.tmpPanel.mask('Loading...');
        try {
            // getting path to ExtJS loading.gif
            SeqView.ExtIconLoading = getComputedStyle(document.querySelector('.x-mask-msg-text'))['background-image'];
        } catch(e) {}
        var options = {url: url, context: this, data: params};
        if (this.m_App.m_parallelRender) options.data.parallel_render = this.m_App.m_parallelRender;

        var checkJobStatus = function(data, text, res) {
            var cleanPanel = function(panel) {
                if (panel.tmpPanel) {
                    delete panel.tmpPanel;
                    panel.removeAll();
                }
            }
            if (this.m_App.m_ReqNum > req_num) return;
            var from_cgi = SeqView.decode(data);
            if (typeof from_cgi !== 'object') {
                cleanPanel(this.m_App.m_Panel);
                fail_fn.apply(scope, arguments);
                return;
            }
            if (from_cgi.job_status) {
                var st = from_cgi.job_status;
                if (st == 'submitted' || st == 'running' || st == 'pending') {
                    options.url = url + '?job_key=' + from_cgi.job_id;
                    Ext.defer(SeqView.App.simpleAjaxRequest,
                        Math.min(15000, Math.max(2000, 1000 * delay++)),
                        this, [options]);
                    return;
                }
            }
            cleanPanel(this.m_App.m_Panel);
            this.m_App.watchdogStop();
            if (from_cgi.success === false || from_cgi.job_status == 'failed') {
                if (!this.m_App.m_NoDataCookie && parKey ){
                    var cookie_name = SeqView.Cookies.UserDataCookieName;
                    if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserDataCookieNameBase) {
                        cookie_name += '-' + this.m_App.m_AppName;
                    }
                    SeqView.SessionData.set(cookie_name, parKey);
                }
                fail_fn.apply(scope, arguments);
            } else {
                Ext.apply(from_cgi.Options, window.JSON.parse(localStorage.getItem('NCBI/SV/Preferences')));
                Ext.apply(this, from_cgi);
                var wrn = '';
                Ext.each(this.SeqInfo.warning_messages || [], function() { wrn += this + '<br>'; });
                this.SeqInfo.warnTT = wrn;
                this.SeqInfo.icon = SeqView.ExtIconLoading.slice(5, -18) + (wrn ? 'shared/warning.gif' : 'window/toast/icon16_info.png');
                succ_fn.call(scope, this);
            }
        }

        this.m_App.m_timeStamps = {prior2seqconfig: new Date().getTime()};
        if (window.timeStamp) {
            this.m_App.ping({'SV_start_to_seqconfig_time': this.m_App.m_timeStamps.prior2seqconfig - timeStamp,
                          'sv-event':'initialization'});
            delete timeStamp;
        }
        
        options.success = options.error = checkJobStatus;
        this.m_App.watchdogStart(url, '', params);
        this.m_App.AjaxRequest(options);
    },

    visualOptionsUrl: function() {
        var params = {};
        if (this.Options) {
            Ext.each(['label','color', 'decor', 'spacing', 'content'], function(p) {
                var pattr = 'curr_'+p;
                if (this.Options[pattr] && this.Options[pattr].length > 0)
                    params[p] = this.Options[pattr];
            },this);
        }
        return params;
    },
    
    getObjIdTrack: function(id) {
        if (!id) return null;
        var str = '';
        Ext.each(this.TrackConfig,
            function(track, ix) {
                if (track.id != id) return true;
                str = 'id:' + track.id + ',' + 'key:' + track.key;
                if (track.name && track.name.length > 0)
                    str += ',' + 'name:' + NCBIGBUtils.escapeTrackName(track.name);
                if (track.data_key && track.data_key.length > 0)
                    str += ',' + 'data_key:' + NCBIGBUtils.escapeTrackName(track.data_key);
                if (track.display_name && track.display_name.length > 0)
                    str += ',' + 'display_name:' + NCBIGBUtils.escapeTrackName(track.display_name);
                if (track.dbname && track.dbname.length > 0 ){
                    str += ',' + 'dbname:' + NCBIGBUtils.escapeTrackName(track.dbname);
                }
                Ext.each(track.annots, function(annot, idx){
                    str += (idx == 0 ? ',annots:' : '|') + NCBIGBUtils.escapeTrackName(annot);
                });

                return false;
            }
        );
        return '[' + str + ']';
    },

    cleanGeneMode: function(t) {
        var vs = t.choice_list[0].values
        var name = null;
        for (var i = vs.length; i--;) {
            if (vs[i].display_name.charAt(0) != '*') continue;
            vs[i].display_name = vs[i].display_name.replace(/^\* /,'');
            name = vs[i].name;
        }
        return name;
    },

    changeGeneMode: function(showAll) {
        Ext.each(this.TrackConfig, function(t) {
            if (t.category.name != 'Genes') return;
            var cur = t.choice_list[0];
            if (showAll) {
                if (Ext.each(cur.values, function(v) { return v.name != 'ShowAll' }) === true) return;
                Ext.each(cur.values, function(v) {
                    if (v.name !== cur.curr_value) return;
                    v.display_name = '* ' + v.display_name;
                    return false;
                });
                cur.curr_value = 'ShowAll';
            } else cur.curr_value = this.cleanGeneMode(t) || cur.curr_value;
        }, this);
    },

    removeUserTrack: function(uuid) {
        var tIdx = -1;
        Ext.each(this.TrackConfig, function(t, idx) {
            if (t.uuid != uuid) return true;
            tIdx = idx;
            return false;
        });
        if (tIdx < 0) return;
        TMS.GetTracksById([this.TrackConfig[tIdx].id]).done(function(){
            var trx = this.GetTracks();
            if (trx.length == 0) tIdx = -1;
            else Ext.each(trx, function() { this.Deregister(); });
        });
        if (tIdx < 0) return;
        this.TrackConfig.splice(tIdx, 1);
    }
};

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Api
/********************************************************************/
SeqView.Api = (function() {

    function constructor(app) {
        this.m_App = app;
    }

    return constructor;

}) ();

SeqView.Api.prototype = {
    getAllDisplayOptions: function() {
        return this.m_App.m_Config.DisplayOptions;
    },

    getTracks: function(all) {
        return this.m_App.getTrackObjects(all);
    },

    setTracks: function(tracks) {
    }
};

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.App
/********************************************************************/

Ext.define('SeqView.App', {
    extend: 'Ext.util.Observable',

    constructor: function(div_id) {
        var sm_ResizeWatcher = null;
        if (this.m_domDiv = Ext.get(div_id))
            this.m_domDiv.addCls("SeqViewerApp");
        else
            throw "A div element containing SeqViewer should have an id attribute";

        var onWindowResize = function(w, h) {
            var rw = w;
            var rh = h;
            if (sm_ResizeWatcher) {
                clearTimeout(sm_ResizeWatcher);
                sm_ResizeWatcher = null;
            }
            sm_ResizeWatcher = Ext.defer(function(w, h) {
                SeqView.m_Apps.forEach(function(app) { app.doWindowResize(w, h); });
            }, 500,this,[rw,rh]);
        }

        Ext.on('resize', onWindowResize);

        this.m_DivId = div_id;
        this.m_DivCtrlId = this.m_DivId + 'Controls';
        this.m_DivJSId = this.m_DivId + 'JS';
        this.m_DivTitle = this.m_DivId + 'Title';
        this.m_DivTitleID = this.m_DivId + 'TitleID';
        this.m_Views = [];
        this.m_Reflections = null;
        this.m_MarkersInfo = null;
        this.m_TextView = null;
        this.m_InitialLoading = false;
        this.m_DialogShown = false;
        this.m_SelectionTTMap = new Map();
        this.m_Idx = SeqView.m_Apps.push(this) - 1;
        this.callParent(arguments);
    },

    pingClick: function(area, svevent) {
            var obj = {jsevent: 'click', 'sv-area': area};
            if (svevent) obj['sv-event'] = svevent;
            this.ping(obj);
    },
    
    tmp_pingArgs: [],
    ping: function(a) { // saving pings while NCBI instrumented page is loading
        var args = this.tmp_pingArgs;
        if (typeof ncbi !== 'undefined' && ncbi.sg.ping) {
            this.ping = function(obj) {
                obj.ncbi_app = 'sviewer-js';
                obj['sv-appname'] = this.m_AppName;
                ncbi.sg.ping(obj);
            };
            if (SeqView.doInitPing) {
                this.ping({'sv-browser': Ext.browser.name,
                    'sv-browser-version': Ext.browser.version.major, 'sv-os': Ext.os.name});
                delete SeqView.doInitPing;
            }
            while (args.length) this.ping(args.shift());
            delete this.tmp_pingArgs;
            this.ping(a);
        } else {
            if (args.push(a) > 5) console.log('DROPPED PING: ' + JSON.stringify(args.shift()));
        }
    },
    
    getApi: function() { return new SeqView.Api(this); },
    getWidth: function() { return this.m_domDiv.getWidth() - (this.m_Embedded ?  0 : 2); },

    showPrintPageDlg: function(app_idx) {
        var app = this.findAppByIndex(app_idx);
        if (app) { app.showPrintPageDlg(); }
    },

    doWindowResize: function(w,h) {
        if (this.m_BrowWidth != w) {
            this.reloadAllViews();
            this.m_BrowWidth = w;
        }
        this.updateLayout();
    },

    reload: function(rel_attr, keep_config) {
        if (this.m_InitialLoading){
            if (this.m_deferredReload) clearTimeout(this.m_deferredReload);
            this.m_deferredReload = Ext.defer(this.reload, 100, this, [rel_attr, keep_config]);
            return;
        }
        delete this.m_deferredReload;
        if (this.m_AjaxTrans) {
            // kill all active ajax requests
/*            try { //doesn't work while we use jQuery AJAX
                Ext.each(this.m_AjaxTrans,function(tid) {
                    Ext.Ajax.abort(tid);
                });
            } catch(e) {}*/
            delete this.m_AjaxTrans;
        }
        this.m_Panorama = null;
        this.forEachView(function(v) { v.destroy(); });
        this.m_Views = [];

        var div = Ext.get(this.m_DivJSId);
        if (div) div.remove();
        div = Ext.get(this.m_DivCtrlId);
        if (div) div.remove();

        if (keep_config !== true)
            delete this.m_Config;
        if (this.m_PermConfId) {
            var el = Ext.get(this.m_PermConfId);
            if (el) {
                var fc = el.first()
                if (fc)
                    fc.remove();
            }
        }
        this.load(rel_attr);
    },

//////////////////////////////////////////////////////////////////////////
// initLinks: initialize GI and data independent links, all other links
// are initialized in loadAccession
    initLinks: function() {
        var prefix = SeqView.base_url;
        var html_prefix = prefix;
        if (html_prefix.length > 0 && html_prefix.charAt(html_prefix.length-1) !== '/') html_prefix += '/';
        this.m_CGIs = {
            prefix:     prefix,
            html_prefix: html_prefix,
            FeatSearch: prefix + 'sv_search.cgi',
            ObjInfo:    prefix + 'seqgraphic.cgi',
            NetCache:   prefix + 'ncfetch.cgi',
            Feedback:   prefix + 'feedback.cgi',
            Config:     prefix + 'seqconfig.cgi',
            Panorama:   prefix + 'seqgraphic.cgi',
            Graphic:    prefix + 'seqgraphic.cgi',
            Alignment:  prefix + '../msaviewer/alnmulti.cgi',
            ObjCoords:  prefix + 'seqgraphic.cgi',
            Sequence:   prefix + 'seqfeat.cgi',
            SequenceSave: prefix + 'sequence.cgi',
            Watchdog:   prefix + 'sv_watchdog.cgi',
            SearchMru:  prefix + 'sv_mru.cgi'
        };
        if (this.m_hiddenOptions.search('use_jsonp') >= 0)SeqView.jsonType = 'JSONP';
        else {
            SeqView.jsonType = 'JSON';
            if (Ext.isIE) {
                SeqView.App.simpleAjaxRequest({url: this.m_CGIs.Feedback,
                    success: function(data, text) { SeqView.jsonType = 'JSON'; },
                    error: function(data, text) { SeqView.jsonType = 'JSONP'; }
                });
            }
        }
        
        if (!this.m_AppName || this.m_AppName == 'no_appname') {
            if (!SeqView.origHostName && !this.m_Embedded && SeqView.standalone) this.m_AppName = 'ncbi_sviewer';
            else {
                var newAppName = document.location.hostname;
                if (SeqView.origHostName || (newAppName.search(/www|qa|test|blast/) != 0)) {
                    var msg = 'Sequence Viewer was initialized with parameter <b>appname=' + this.m_AppName + '</b><br>';
                    Ext.Msg.show({title: 'Warning',
                        msg: msg + 'It will be substituted with the hostname: <b>' + newAppName + '</b>',
                        icon: Ext.Msg.WARNING,
                        buttons: Ext.Msg.OK
                    });
                }
                this.m_AppName = newAppName;
                this.m_appWarning = 'Uninitialized <b>appname</b> is set to <b>' + newAppName + '</b><br>';
           }
        }
     },

    load: function(rel_attr, maxDelay) {
        if (this.m_Config) {
            console.warn('SeqView.App.load() can be invoked only once per app instance. Use SeqView.App.reload()');
            return;
        }
        maxDelay = maxDelay || 1000;
        var interval = 3;
        var app = this;
        setTimeout(function() { // Delay to allow jQuery to get ready
            maxDelay -= interval;
            if (typeof TMS === 'undefined' && maxDelay > 0)
                app.load(rel_attr, maxDelay);
            else {
                app.parseParams(rel_attr);
                app.initLinks();
                app.m_Config = new SeqView.Config(app);
                app.loadExtraData();
            }
        },interval);
    },

    // callback object with success and failure defined
    reset: function(callback, reload) {
        this.m_Config = new SeqView.Config(this);
        if (reload) {
            if (this.m_MarkersInfo) {
                this.m_MarkersInfo.reset();
            }
            if (this.m_Reflections) {
                this.m_Reflections.deleteAll();
            }
            this.forEachView( function(view) { view.remove(); } )
            this.m_Views = [];
            this.m_InitialLoading = true;
        }
        var succ_fn = callback && callback.success ? callback.success : function() {};
        var fail_fn = callback && callback.failure ? callback.failure : function() {};
        this.m_Config.load({
            defaultcfg: true,
            resettracks: true,
            callback: {
                success: function(res) { if (reload) { this.infoLoaded(res); } succ_fn(); },
                failure: function(data, text, res) { this.infoMissed(data, text, res); fail_fn(); },
                scope:   this
            }
        });
    },

    // repackage DisplayOptions and TrackConfig to the shape ready for consumption
    // by api functions getDisplayOptions and getTracks
    postProcessConfig: function() {
        // Recode m_Config.DisplayOptions from array to map from track key to track object
        if (!('DisplayOptions' in this.m_Config)) {
            return;
        }
        var tracks = this.m_Config.DisplayOptions.tracks;
        var new_tracks = {};
        for (var i = 0; i < tracks.length; i++) {
            var key = tracks[i].key;
            if ('subkey' in tracks[i]) {
                key += '.' + tracks[i].subkey;
            }
            new_tracks[key] =  tracks[i];
        }
        this.m_Config.DisplayOptions.tracks = new_tracks;
    },

    getTrackObjects: function(all) {
        var tracks = [];
        Ext.each(this.m_Config.TrackConfig, function(track) {
            if (track.shown || all) {
                tracks.push(track);
            }
        });
        tracks.sort(function(t1,t2) {
            return t1.order - t2.order;
        });
        return tracks;
    },

    getTracks: function() {
        return this.getActTracks();
    },

    getActTracks: function() {
        if (this.m_actTracks == '') {
            this.m_actTracks = SeqView.TM.tracksArrayToString(this.getActiveTracks(), true);
        }
        return this.m_actTracks;
    },

    getConfiguration: function() {
        return this.m_Config;
    },

    // Returns an array of track info objects intended for easy constructing
    // 'tracks' URL parameter for subsequence embedding of Sequence Viewer
    // For complete track string either call getTrackStringFromInfo, or embed
    // it's code into your own app.
    getTracksInfo: function() {
        tracks = [];
        for (var i = 0; i < this.m_Config.TrackConfig.length; i++) {
            var track = this.m_Config.TrackConfig[i];
            tracks.push({
                name:  track.name,
                order: track.order,
                shown: track.shown,
                descr: SeqView.TM.trackToString(track)
            });
        }
        return tracks;
    },

    // Returns tracks parameter for Sequence Viewer, use it to construct parameter
    // string for embedding. Does not have parameter name, 'tracks', only value -
    // series of brace enclosed strings describing tracks.
    getTrackStringFromInfo: function(tracks_array) {
        var str = '';
        for (var i = 0; i < tracks_array.length; i++) {
            var track = tracks_array[i];
            str += '[' + track.descr;
            str += ',order:' + track.order;
            str += ',shown:' + (track.shown ? 'true' : 'false');
            str += ']';
        }
        return str;
    },

//////////////////////////////////////////////////////////////////////////
// clearParams:

    clearParams: function() {

        this.GI = null;
        this.m_From = null;
        this.m_To = null;
        this.m_ViewSearch = null;
        this.m_ViewSlimModes = [];
        // New markers
        this.m_Markers = [];
        this.m_Views = [];
        this.m_Reflections = null;
        this.m_MarkersInfo = null;
        this.m_TextView = null;

        this.m_ViewMarkers = false;
        this.m_Key = null;
        this.m_NAA = null;
        this.m_BamPath = null;
        this.m_SRZ = null;
        this.m_hiddenOptions = '';
        this.m_actTracks = "";
        this.m_userData = [];
        this.m_userTrackNames = [];
        this.m_FindCompartments = true;
        this.m_NeedAlignmentView = false;
        this.m_OnlyAlignmentView = false;
        this.m_DepthLimit = null;
        this.m_AppContext = null;

        this.m_Origin = 0;
        this.m_Flip = false;
        this.m_ViewLabels = null;
        this.m_Embedded = false;
        this.m_SnpFilter = null;
        this.m_ViewSeq = null;

        this.m_ViewOptions = {};
        this.m_ViewRanges = [];
        this.m_ViewColors = [];
        this.m_ViewsSelect = [];
        this.defaultConfig = {};
        delete this.m_ViewTheme;
        delete this.m_ViewContent;
        delete this.m_AppName;

        delete this.m_Panel;
        this.m_Panel = null;
        this.m_ViewParams = null;
        this.m_ItemID = null;

        this.m_LocalPrimerTool = false;
        this.m_LocalBlast = false;
        this.m_ShowQuery = false;
        this.m_QueryRange = null;
        this.m_NoDataCookie = false;

        delete this.m_highlightsColor;
        this.m_GraphicExtraParams = {};

        this.m_Toolbar = {
            history: true,
            name: true,
            search: true,
            panning: true,
            zoom: true,
            genemode: true,
            slimmode: true,
            tools: true,
            config: true,
            reload: true,
            help: true
        };
    },

    parseParams: function(rel_attr) {
        this.clearParams();
        delete SeqView.requestTrackSets;
        var doc_location = document.location.href;
        this.m_AllViewParams = doc_location.split("?")[1] || '';
        this.m_AllViewParams = this.m_AllViewParams.replace(/#/, ''); // safety

        if (!rel_attr) {
            var rel = Ext.get(this.m_DivId);
            if (rel && rel.dom) {
                var aref = rel.first();
                if (aref && aref.dom && aref.dom.tagName == "A") {
                    rel_attr = aref.dom.href.slice(aref.dom.href.indexOf('?') + 1)
                    aref.remove();
                } else { // for compatibility with the old verions
                    rel_attr = rel.dom.getAttribute('rel');
                }
            }
        }
        if (rel_attr != null) {
            var href = this.m_AllViewParams;
            var keep = ''; 
            var pLst = ['parallel_render', 'extra_opts'];
            Ext.each(pLst, function(p){
                var idx = href.indexOf(p);
                if (idx < 0) return;
                keep += '&' + href.slice(idx).split('&')[0]
            });
            // We should not merge URL parameters in for embedded viewers.
            // Better is to parse rel_attr into parameter set first,
            // but this is reliable enough - equal sign in all other
            // places beyond parameter definition should be escaped
            this.m_AllViewParams = (/embedded=/.test(rel_attr))
                ? rel_attr
                : this.m_AllViewParams + '&' + rel_attr;
            this.m_AllViewParams += keep;
        }
        this.m_AllViewParams = Ext.util.Format.htmlDecode(this.m_AllViewParams);

        this.m_MarkersInfo = new SeqView.MarkersInfo(this);

        var marker_ranges = [];
        var marker_names = [];
        if (this.m_AllViewParams) {
            var p_array = this.m_AllViewParams.split('&');

            delete this.mf_MultiPanel;
            var sKeys= ['key', 'naa', 'tracks', 'tkey', 'tn', 'rkey', 'mn', 'mk',
                        'url', 'url_reload', 'rid', 'data', 'id', 'select',
                        'srz', 'content', 'viewer_context', 'app_context',
                        'assm_context', 'snp_filter', 'bam_path', 'theme', 'search'];

            for (var i = 0; i < p_array.length; i++) {
                var pair = p_array[i]
                var p = pair.split("=");
                var the_key = p[0].toLowerCase();

                var val = (p[1] === undefined ? '' : unescape(p[1]));
                
                if (sKeys.indexOf(the_key) == -1) val = val.toLowerCase();
                
                switch (the_key) {
                    case 'bam_path': this.m_BamPath = val; break;
                    case 'depth_limit': this.m_DepthLimit = val;  break;
                    case 'id': this.GI = val; break;
                    case 'f': case 'from': this.m_From = val; break;
                    case 't': case 'to':  this.m_To = val; break;
                    case 'gl_debug': if (val == '1' || val == 'true') this.m_GraphicExtraParams['gl_debug'] = 'true'; break;
                    case '_gene': case 'gene': this.m_ViewSearch = val; break;
                    // old way to parse markers
                    case 'm': marker_ranges = val.split(','); break; // like: m=249,276
                    case 'mn': marker_names = val.split(','); break; // like: mn=Marker 1,Marker 2
                    // new way, integrated parameter mk=249:250|Marker 1|ff0000,276|Marker 2|00ff00
                    case 'mk': this.m_MarkersInfo.parseMarkersURL(this.m_Markers, val);  break;
                    case 'key': this.m_Key = val; break;// like: db=nucleotide
                    case 'naa': this.m_NAA = val; break; // like: naa=NA000000003,NA000000004, NA000000004:renderer_name
                    case 'search': this.m_ViewSearch = val; break; // like: search=feature_name
                    case 'srz': this.m_SRZ = val; break; // like: srz=SRZ000200
                    case 'tracks': this.defaultConfig.m_TracksFromURL = this.m_TracksFromURL = val; break;
                    case 'nocache': this.defaultConfig.nocache = val; break;
                    case 'vm': this.m_ViewMarkers = true; break; // like: vm=true
                    case 'tn': this.m_userTrackNames.push(decodeURIComponent(val)); break;
                    case 'find_comp': this.m_FindCompartments = val == "true" || val == "1" || val == "on"; break;
                    case 'data': val = decodeURIComponent(val); // literal data, URI component encoded
                    case 'rkey':
                    case 'rid':  
                    case 'url_reload': // the same like 'url' with check_cs = true 
                    case 'url': this.m_userData.push([the_key, val]); break; // like: url=www.ncbi.nlm.nih.gov/data.txt
                    case 'align': this.m_NeedAlignmentView = true; break;
                    case 'onlyalign': this.m_OnlyAlignmentView = true; break;
                    case 'labels': this.m_ViewLabels = val; break;
                    case 'embedded': this.m_Embedded = val == "true" || val == "1" || val; break;
                    case 'multipanel': this.mf_MultiPanel = (val == '1' || val == 'true' || val == '' || val == null); break;
                    case 'snp_filter': this.m_SnpFilter = val; break;// like: Validated|000010000%20-1%20-1%20-1%20_%20_%20_%2043%20_%20_%20
                    case 'origin': this.m_Origin = parseInt(val); break;// like: origin=10000 ; 0-based
                    case 'seq': this.m_ViewSeq = val.split(':'); break;// like: seq=2000:7000

                    case 'highlights_color': this.m_highlightsColor = val; break;// default color for highlights
                    case 'color': // Options are global to SeqView app take only 0-item
                    case 'label': 
                    case 'decor': 
                    case 'spacing': this.m_ViewOptions[the_key] = val.split(',')[0]; break;
                    case 'theme': this.m_ViewTheme = val.split(',')[0]; break; // like: theme=Compact,Details
                    case 'content': this.m_ViewContent = val.split(',')[0]; break;

                    case 'v': if(val.length) { this.m_ViewRanges = val.split(','); } break; // like: v=1000:6000,8000:18723
                    case 'c': this.m_ViewColors = val.split(','); break; // like: color=0000FF,FFFF99
                    case 'slim': this.m_ViewSlimModes = val.split(','); break; // like: slim=1,0,... or slim=true,false,...
                    case 'gflip': this.m_Flip = val=='true' || val=='1'; break; // gflip=true
                    case 'flip': case 'strand': var parts = val.split(','); // like: flip=true,false
                        if (parts && parts.length > 0) this.m_Flip = parts[0]=="true" || parts[0]=="1";
                        break;
                    case 'select': this.m_ViewsSelect = val.split(','); break; // like: select=
                    case 'itemid': this.m_ItemID = val; break;
                    case 'noviewheader': this.m_NoViewHeader = val=='true' || val=='1'; break;
                    case 'viewer_context': this.m_ViewerContext=val; break;
                    case 'app_context': this.m_AppContext=val; break;
                    case 'assm_context': this.m_AssmContext=val; break;
                    case 'noguess_assm': this.m_noGuessAssm = val=='true' || val=='1'; break;
                    case 'nopdf': this.m_NoPDF = val=='true' || val=='1'; break;

                    case 'extra_opts': this.m_hiddenOptions = val; break;
                    case 'parallel_render': this.m_parallelRender = val; break;
                    case 'iframe': this.m_iFrame = val; break;
                    case 'nodatacookie': this.m_NoDataCookie = val=='true' || val=='1'; break;
                    case 'showquery': this.m_ShowQuery = val=='true' || val=='1'; break;
                    case 'queryrange': var parts = val.split(':');
                        if (parts && parts.length > 1) {
                            this.m_ShowQuery = true;
                            this.m_QueryRange = [ parseInt(parts[0]), parseInt(parts[1])];
                        }
                        break;

                    case 'primer': if (location.pathname.indexOf('staff') != -1) this.m_LocalPrimerTool = val=='true' || val=='1'; break;
                    case 'blast': if (val == 'web') this.m_LocalBlast = 'web'; else this.m_LocalBlast = val=='true' || val=='1'; break;
                    case 'tkey': this.m_TKey = val; break;
                    case 'appname': if (val && val.length > 0) this.m_AppName = val; break;
                    case 'noconfdlg': this.m_NoConfDlg = val=='true' || val=='1'; break;
                    case 'nosavetracks': this.m_noSaveTracks = true; break;
                    case 'toolbar': if (!val) break;
                        var add = val.charAt(0) != '-';
                        var start = add ? 0 : 1;
                        if (add) this.m_Toolbar = {};
                        for (var j = start; j < val.length; j++) {
                            switch (val.charAt(j)) {
                                case 'b': this.m_Toolbar["history"] = add; break;
                                case 'n': this.m_Toolbar["name"] = add; break;
                                case 's': this.m_Toolbar["search"] = add; break;
                                case 'p': this.m_Toolbar["panning"] = add; break;
                                case 'z': this.m_Toolbar["zoom"] = add; break;
                                case 'g': this.m_Toolbar["genemode"] = add; break;
                                case 'm': this.m_Toolbar["slimmode"] = add; break;
                                case 't': this.m_Toolbar["tools"] = add; break;
                                case 'c': this.m_Toolbar["config"] = add; break;
                                case 'r': this.m_Toolbar["reload"] = add; break;
                                case 'h': this.m_Toolbar["help"] = add; break;
                            }
                        }
                        break;
                }
                SeqView.TM.renderStat = this.m_hiddenOptions.indexOf('render_stat') >= 0;
            } // for

            if( this.mf_MultiPanel === undefined ){
                this.mf_MultiPanel = !this.m_Embedded;
            }

            if (this.m_Embedded && !this.m_NoDataCookie) {
                var cookie_name = SeqView.Cookies.UserDataCookieName;
                if (this.m_AppName && this.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserDataCookieNameBase) {
                    cookie_name += '-' + this.m_AppName;
                }

                if( !this.m_Key ){
                    // if the NC key with the user's data is not passed check if it is set in cookies
                    this.m_Key = SeqView.SessionData.get( cookie_name, null );
                } else {
                    // save it into a session cookie
                    SeqView.SessionData.set( cookie_name, this.m_Key );
                }
            }
        } // view
        this.m_MarkersInfo.parseMarkersOldStyle(this.m_Markers, marker_ranges, marker_names);

        this.m_Portal =
            this.m_Embedded === false
            && (this.m_AppName || '').indexOf('ncbientrez') == 0;

        var tmpl = new Ext.Template(
            '<a onClick="SeqView.App.showLinkURLDlg(' + this.m_Idx + ', \'3-0\');" href="#">Link To This View</a> | ',
            '<a onClick="SeqView.App.showFeedbackDlg(' + this.m_Idx + ', \'3-1\');" href="#">Feedback</a>');
        var cfg = [
            {tag: 'div', cls: 'SeqViewerControls hidden_for_print', id: this.m_DivCtrlId, html:tmpl.apply()},
            {tag: 'div', cls: 'SeqViewerJS', id: this.m_DivJSId}];

        if (this.m_Embedded) {
            if (!Ext.get(this.m_PermConfId = this.m_DivId + '_confdlg'))
                delete this.m_PermConfId;
            cfg.shift();
            cfg[0].cls = '';
        }
        if( this.GI || this.m_userData.length ){
            Ext.DomHelper.append(this.m_DivId, cfg);
            var view_config = {
                collapsible: (this.m_Embedded === false),
                renderTo: this.m_DivJSId,
                html: '<span id="string_ruler_unit"></span>',
                header: false
            };
            this.m_Panel = new Ext.Panel(view_config);
            this.m_Panel.tmpPanel = this.m_Panel.add(new Ext.Panel({
                height: 60, border: false,
                items:[{xtype:'displayfield', name: 'progress', id:'sv-uplmsg' + this.m_Idx, value: '',
                             style: {color:'grey', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}},
                       {xtype:'displayfield', name: 'error', id:'sv-uplerr' + this.m_Idx, value: '',
                             style: {color:'red', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}}]
            }));
            this.m_Panel.updateLayout();
        }
    },

    resizeIFrame: function(height){
        if (!this.m_iFrame) return;
        var h = document.documentElement;
 /*       var b = document.body;
        var svH = this.m_domDiv.getHeight();
        console.log(this.m_iFrame, svH, height || 0, b.scrollHeight, b.offsetHeight, h.clientHeight, h.scrollHeight, h.offsetHeight);*/
        var height = Math.max(height || 0, h.offsetHeight + (Ext.isIE ? 9 : Ext.isFirefox ? 5 : 0));
        window.frames.parent.postMessage({f: this.m_iFrame, h: height}, '*');
    },

//////////////////////////////////////////////////////////////////////////
// openNewWindowPOST(url, params) - open new window by POSTing 'params' to 'url'
//     url - URL to use for POST
//     params - object with parameters to pass in the POST request
    openNewWindowPOST: function(url, params) {
        var new_win = window.open('about:blank', '_blank');
        new_win.focus();
        var new_body = new_win.document['body'];
        var form = new_win.document.createElement('form');
        form.method = 'POST';
        for (var name in params) {
            if (!params.hasOwnProperty(name)) continue;
            var el = new_win.document.createElement('input');
            el.type = 'hidden';
            el.name = name;
            el.value = params[name];
            form.appendChild(el);
        }
        new_body.appendChild(form);
        form.action = url;
        form.submit();
    },

    blast: function(range) {
        var params = {};
        params["SV_SHOW"] = "TRUE";
        if (this.m_Key) {
            params["SV_KEY"] = this.m_Key;
        }
        params["PAGE_TYPE"] = "BlastSearch";
        params["SHOW_DEFAULTS"] = "on";
        if (this.m_ViewParams['acc_type']=='protein') {
            params["PROGRAM"] = "blastp";
            params["BLAST_PROGRAMS"] = "blastp";
        } else {
            params["PAGE"] = "Nucleotides";
            params["PROGRAM"] = "blastn";
            params["MEGABLAST"] = "on";
            params["BLAST_PROGRAMS"] = "megaBlast";
        }
        if (range) {
            params["QUERY_FROM"] = (range[0]+1);
            params["QUERY_TO"]   = (range[1]+1);
        }
        var prod = SeqView.webNCBI.indexOf("www.ncbi.nlm.nih.gov") != -1;
        var prefix = this.m_LocalBlast === 'web' ? "https://web.ncbi.nlm.nih.gov/blast/" : (
                                           this.m_LocalBlast ? "" : (prod ? "https://blast.ncbi.nlm.nih.gov/" :
                                           SeqView.webNCBI + "blast/"));
        var url = prefix + "Blast.cgi";
        // check local id, fill out query with sequence in this case
        if (this.GI.indexOf("lcl|") != 0) {
            params["QUERY"] = this.GI;
            this.openNewWindowPOST(url, params);
        } else {
            var ranges = "0-" + (this.m_SeqLength-1);
            this.AjaxRequest({
                url: this.m_CGIs.SequenceSave, context: this,
                data: { id: this.GI, ranges: ranges, format: "fasta", nofile: "true", key: this.m_Key },
                dataType: "text",
                success: function(data, text, res) {
                    // The define line is coming with coordinate range, which Blast does not understand.
                    // So, we need to remove it from, e.g. ">lcl|Query_1:1-2534" to read ">lcl|Query_1"
                    params["QUERY"] = data.replace(/>(lcl\|.*):\d+-\d+/, function(m, p1) { return ">" + p1; });
                    params["QUERY_BELIEVE_DEFLINE"]="true";
                    this.openNewWindowPOST(url, params);
                },
                error: function(data, text, res) {
                    console.log("sequence.cgi failure - " + res);
                }
            });
        }
    },

//////////////////////////////////////////////////////////////////////////
// primerBlast:
    primerBlast: function(whole, range_set) {
        var seqinfo = this.m_ViewParams;
        if (seqinfo["acc_type"] !== "DNA") return;
        if (!whole  &&
            (!range_set || range_set.length < 1 || range_set.length > 2)) return;
        var url = this.m_LocalPrimerTool
            ? "primer-blast/index.cgi?"
            : "/tools/primer-blast/index.cgi?"
        ;
        url += "ORGANISM=" + encodeURI(seqinfo["organism"]);
        url += "&INPUT_SEQUENCE=" + encodeURI(seqinfo["id"]);
        var total_range;
        if (whole)
            total_range = [0, this.m_SeqLength-1];
        else
            total_range = range_set[0];
        if (!whole  &&  range_set.length === 2) {
            var range0, range1;
            if (range_set[0][1] < range_set[1][0]) {
                var range0 = range_set[0];
                var range1 = range_set[1];
            } else {
                var range0 = range_set[1];
                var range1 = range_set[0];
            }
            url += "&PRIMER5_END=" + (range0[1]+1);
            url += "&PRIMER3_START=" + (range1[0]+1);
            total_range = [ range0[0], range1[1] ];
        }
        url += "&PRIMER5_START=" + (total_range[0]+1);
        url += "&PRIMER3_END=" + (total_range[1]+1);
        url += "&PRIMER_PRODUCT_MAX=" + (total_range[1] - total_range[0]);
        url += "&SHOW_SVIEWER=true";
        if (this.m_Key) {
            url += "&SVIEWER_DATA_KEY=" + this.m_Key;
        }
        window.open(url);
    },

    createAndLoadMainPanel: function() {
        this.m_BrowWidth = Ext.getBody().getViewSize().width;
        Ext.fly(this.m_DivId).on({'contextmenu': this.onContextMenu, scope: this});
        this.loadAccession();
    },

    onContextMenu: function(e) {
        if (e.getTarget().id.indexOf('combobox')) e.stopEvent();
    },

    requestGI: function() {
        Ext.get('app-frontpage').setStyle('display', 'block');
        var button = Ext.get('mainshow-btn');
        button.on('click', function() { this.loadAccessionDlg();  },this);
    },

    loadAccessionDlg: function() {
        Ext.MessageBox.prompt('Load Accession', 'Please enter accession or GI:', function(btn, text) {
            if (btn!='ok'  || text.length==0) return;
            this.redirectWithGi(text);
        }, this, false, this.GI);
    },

    redirectWithGi: function(text_gi) {
        var url = this.m_CGIs.Config;
        this.AjaxRequest({url:url  + '?id=' + (text_gi) + '&notrackconfig=true', context: this,
            success:function(data, text, res) {
                var from_cgi = SeqView.decode(data);
                if (from_cgi.job_status) {
                    var st = from_cgi.job_status;
                    if (st == 'submitted' || st == 'running' || st == 'pending') {
                        Ext.defer(SeqView.App.simpleAjaxRequest, 2000, this, {url: url + '?job_key=' + from_cgi.job_id, context: this});
                        return;
                     }
                }
                if (from_cgi.error) {
                    Ext.MessageBox.show({title: 'Sequence Viewer', msg:from_cgi.error,
                                         buttons: Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                } else if(from_cgi.success === false) {
                    Ext.MessageBox.show({title: 'Sequence Viewer', msg:from_cgi.msg,
                                         buttons: Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                } else {
                    var href = document.location.href.replace(/#/, '');
                    var pars = href.split("?");
                    var args = pars[1] ? pars[1].split('&') : [];
                    var view_url = pars[0] + '?id=' + text_gi;
                    Ext.each(args, function(v) {
                        if (v.indexOf('id=') != 0) view_url += '&' + v;
                    });
                    document.location.href = view_url;
                }
            }});
    },

    showTracksConfigDlg: function(area) {
        this.pingClick(area);
        SeqView.TM.ShowConfigDialog(this, 0);
    },

    updateConfig: function(newconfig, user_data, tracks) {
        this.m_Config.TrackConfig = newconfig.TrackConfig;
        Ext.apply(this.m_Config.Options, newconfig.Options);
        this.m_actTracks = tracks;
        this.reloadAllViews();
        delete this.m_currentTrackSetId;
        this.fireEvent('configuration_changed', this);
    },


//////////////////////////////////////////////////////////////////////////
// showTracks: show/hide tracks matching the query. NB: This function does NOT
//             refresh the view, do it yourself after the call if it returned
//             true.
//    query - object to match with the track to operate on, currently 2 keys
//            supported - "category" and "name". Conditions are joned by AND
//            so if you turn on track with {category: "Variation", name: "SNP" }
//            only tracks with matching category AND names become active.
//            Default (empty object) means ALL tracks.
//    show  - boolean, whether to turn tracks on or off, default off
//    save  - should we save new track configuration, default no
//  returns whether track configuration changed
    showTracks: function(query, show, save) {
        var active_tracks = [];
        var changed = false;
        Ext.each(this.m_Config.TrackConfig, function(track) {
            var shown = track.shown;
            // concatenate conditions with AND
            var select = false;
            for (var p in query) {
                var pp = {}; 
                for (pp in query[p]) {
                    if (pp === '0') {
                        select |= track[p] === query[p];
                        break;
                    }
                    select |= track[p][pp] === query[p][pp];
                }
             }
            if (select) track.shown = show;
            if (track.shown) active_tracks.push(track);
            if (shown != track.shown) changed = true;
        });
            
        active_tracks.sort(function(t1,t2) { return t1.order - t2.order; });
        this.m_actTracks = SeqView.TM.tracksArrayToString(active_tracks, true);
        if (save) {
            var tracks_for_save = SeqView.TM.tracksArrayToString(this.m_Config.TrackConfig, true, true);
            // We don't handle callback succes/error here because we already modified the original
            // config, and disregarding whether it saved successfully or not we can't undo this.
            this.m_Config.save({tracks: tracks_for_save});
        }
        return changed;
   },

    getActiveTracks: function(addHidden) {
        var active_tracks = [];
        var subtracks = [];
        var trackConfig = this.m_Config.TrackConfig;
        Ext.each(trackConfig, function(track) {
            if (!track.shown) return true;
            active_tracks.push(track);
            Ext.each(track.subTracks, function(uId, idx) {
                var legend = track.legend[idx];
                if (addHidden || uId.substr(-7) !== '_hidden') subtracks.push(legend);
                Ext.each(trackConfig, function(ct, cfgIdx) {
                    if (ct.id == legend.id) {
                        var opacity = Math.round(legend.color.split(' ')[3]/2.55);
                        Ext.each(ct.hidden_settings, function() { legend[this.name] = this.value; });
                        if (!legend.hasOwnProperty('opacity')) {
                           ct.hidden_settings.push({name: 'opacity', value: opacity});
                           legend.opacity = opacity;
                        }
                        legend.idx = cfgIdx;
                        return false;
                    }
                });
            });
        });
        active_tracks.sort(function(t1,t2) { return t1.order - t2.order; });

        Ext.each(subtracks, function() {
            var trk = trackConfig[this.idx];
            if (trk && !trk.shown) active_tracks.push(trk);
        });
        return active_tracks;
    },

    loadAccession: function() {
        if (!this.GI) return;

        if (this.m_MarkersInfo) this.m_MarkersInfo.reset();
        if (this.m_Reflections) this.m_Reflections.deleteAll();

        this.forEachView(function(view) { view.remove(); })
        this.m_Views = [];
        this.m_InitialLoading = true;

        var cfg = {callback: {success: this.infoLoaded, failure: this.infoMissed, scope: this}};

        if (this.m_TKey && this.m_TKey.length > 0) {
            var params = {what:'tracks', key: this.m_TKey};
            this.AjaxRequest({
                url: UUD.GetUUDCgiUrl('dl-meta'), data: params, context:this,
                success: function(data, text, res){
                    if (data.success === true)
                        this.defaultConfig.m_TracksFromURL = this.m_TracksFromURL = data.tracks;
                    delete this.m_TKey;
                    this.m_Config.load(cfg);
                },
                error: function(data, text, res){
                    delete this.m_TKey;
                    this.m_Config.load(cfg);
                }
            });
        } else {
            this.m_Config.load(cfg);
        }
    },


    infoMissed: function(data, text, res) {
        var frontpage = Ext.get('app-frontpage');
        if (frontpage) frontpage.setStyle('display', 'block');
        var button = Ext.get('mainshow-btn');
        if (button) button.on('click', function() { this.loadAccessionDlg(); }, this);

        this.m_InitialLoading = false;
        if (data === undefined) return;
        var dt = SeqView.decode(data || {});
        this.m_Panel.add({html: 'An internal error has occurred that prevents Sequence Viewer from displaying.<br>'
            + 'Technical details (seqconfig error): '
            + (dt.msg || dt.statusText || dt.error_message || text)});
    },

    infoLoaded: function(seqconfig) {
        var frontpage = Ext.get('app-frontpage');
        if (frontpage) frontpage.remove();

        if (this.m_PermConfId)
            SeqView.TM.ShowConfigPanel(this, 0);
            
        this.m_ViewParams = this.m_Config.SeqInfo;
        this.m_AssmContext = this.m_ViewParams.assm_context;
        this.m_actTracks = SeqView.TM.generateTracksString(SeqView.TM.processTracksInfo(this.m_Config.TrackConfig));

        this.postProcessConfig();

        this.fireEvent('configuration_loaded', this);

        if (this.m_Embedded === false)
            document.title = this.m_ViewParams['id'] + ': ' + this.m_ViewParams['title'];

        var title_place = Ext.get(this.m_DivTitle);
        if (title_place) title_place.update(this.m_ViewParams['title']);

        var title_id = Ext.get(this.m_DivTitleID);
        if (title_id) title_id.update(this.m_ViewParams['id_full']);

        this.m_Panel.setTitle(this.m_ViewParams['id']);
        this.m_SeqLength = this.m_ViewParams['length'];
        var views = this.m_ViewParams['views'];
        var suggested_ranges = [];
        for (var i = 0; i != views.length; i++) {
            var sp = views[i].split(':');
            // To prevent showing only alignment for PrimerBlast we make an exception here
            if (sp[0]=='multialign' && this.defaultConfig.m_ViewContent != 'PrimerBlast') {
                this.m_NeedAlignmentView = true;
            }
            if (sp[0]=='graphical'  &&  sp[1]) {
                suggested_ranges = sp[1].split('-');
            }
            views[i] = sp[0];
        }

        Ext.each(['label','color','decor','spacing'], function(par) {
            var pidx = this.m_ViewOptions[par] || -1;
            if (pidx >= 0 && pidx < this.m_Config.Options.controls[par].length)
                this.m_Config.Options['curr_' + par] = this.m_Config.Options.controls[par][pidx];
        },this);
        
        if( this.m_Embedded === false || this.m_Embedded == 'panorama' || this.m_Embedded == 'full') {

            this.createPanorama();
        }

        var from, to;
        if (this.m_From && (this.m_From=="begin"  ||  this.m_From=="begining")) this.m_From = null;
        if (this.m_To && this.m_To=="end") this.m_To = null;
        if (!this.m_From && !this.m_To) {
            if (suggested_ranges[0]) {
                from =  NCBIGBUtils.stringToNum(suggested_ranges[0]);
                to  =   NCBIGBUtils.stringToNum(suggested_ranges[1]);
                m_From = (from + 1) + ''; // convert to string
                m_To   = (to + 1) + ''; // convert to string
            } else {
                from = 0;
                to = this.m_SeqLength-1;
                // to = this.m_SeqLength < 10000 ? this.m_SeqLength-1 : this.m_SeqLength > 500000 ? 75000 : this.m_SeqLength / 10;
            }
        } else {
            from = this.m_From ? NCBIGBUtils.stringToNum(this.m_From)-1 : 0;
            to   = this.m_To   ? NCBIGBUtils.stringToNum(this.m_To)-1   : this.m_SeqLength-1;
            // to   = this.m_To   ? NCBIGBUtils.stringToNum(this.m_To)-1   :
            //    (this.m_SeqLength < 10000 ? this.m_SeqLength-1 : this.m_SeqLength > 500000 ? 75000 : this.m_SeqLength / 10);
        }

        if (this.m_QueryRange) {
            var overhang = Math.round(0.15 * (this.m_QueryRange[1] - this.m_QueryRange[0]));
            this.m_ViewRanges = [ Math.max(0, this.m_QueryRange[0]-overhang ) + ":" +
                                  Math.min(this.m_SeqLength-1, this.m_QueryRange[1]+overhang)];
        }

        // fix for old Portal format
        if (this.m_ViewRanges.length==1 && this.m_ViewRanges[0].indexOf('begin') != -1){
            this.m_ViewRanges = []; // ignore cases such as begin..end
        }
        if (this.m_Embedded == 'panorama') return;
        if (this.m_ViewRanges.length == 0) {
            this.m_ViewRanges[0] = (from + 1) + ':' + (to + 1);
        }

        if ((this.m_OnlyAlignmentView === true && this.m_Embedded === false) || (this.m_Embedded && this.m_NeedAlignmentView === true)) {
            var alview = new SeqView.Alignment(this);
            var the_r = this.m_ViewRanges[0];
            var range;

            if (the_r.indexOf(':') != -1) range = the_r.split(':');
            if (the_r.indexOf('..') != -1) range = the_r.split('..');
            if (the_r.indexOf('-') != -1) range = the_r.split('-');
            var zs = range[1] == 'zs';
            range = this.decodeRange(range);
            alview.m_FromSeq = range[0];
            //alview.m_LenSeq  = range[1] - range[0] + 1;
            alview.m_LenSeq  = this.m_SeqLength;
            if (zs) {
                alview.m_LenSeq -=  Math.floor(415 * SeqView.MinBpp); // adjust to align viewer borders
            }
            this.registerView(alview);
        } else {
            for (var i = 0; i < this.m_ViewRanges.length; ++i) {
                var view = new SeqView.Graphic(this);
                var the_r = this.m_ViewRanges[i];
                var range;

                if (the_r.indexOf(':') != -1) range = the_r.split(':');
                if (the_r.indexOf('..') != -1) range = the_r.split('..');
                if (the_r.indexOf('-') != -1) range = the_r.split('-');
                var zs = range[1] == 'zs';

                range = this.decodeRange(range);

                view.m_FromSeq = range[0];
                view.m_LenSeq  = range[1] - range[0] + 1;

                view.m_UrlFrom = range[0] ? range[0] +1 : null;
                view.m_UrlTo   = range[1] ? range[1] +1 : null;
                view.m_slimMode = false;
                var vsm = this.m_ViewSlimModes;
                if (vsm.length > i) view.m_slimMode = (vsm[i] == "true" || vsm[i] == "1");
                if (this.m_ViewColors.length > i) view.m_Color = this.m_ViewColors[i];
                view.setFlipNoReload(this.m_Flip, true);

                if (this.m_ViewsSelect.length > i  && this.m_ViewsSelect[i] != 'null') view.m_SelectedSig = this.m_ViewsSelect[i];

                this.registerView(view);

                if (i == 0 && this.m_NeedAlignmentView && this.m_Embedded === false) {
                    var alview = new SeqView.Alignment(this);
                    alview.m_FromSeq = range[0];
                    alview.m_LenSeq  = range[1] - range[0] + 1;
                    if (zs) {
                        alview.m_LenSeq -=  Math.floor(415 * SeqView.MinBpp); // adjust to align viewer borders
                    }
                    this.registerView(alview);
                }
            }
        }

        if (this.m_ViewSeq && this.m_ViewSeq.length > 0) {
            var from = NCBIGBUtils.stringToNum(this.m_ViewSeq[0])-1;
            var to = NCBIGBUtils.stringToNum(this.m_ViewSeq[1])-1;
            this.createTextView([from,to,this.m_Flip]);
        }
    },

    getToolTipForTrack: function(view,tname) {
        //for tracks extracting first part, delimited by ";"
        if(tname.match(";")){
            var arr = [];
            arr = tname.split(";");
            var gijunk = arr[1];//arr[1] is used for getting track name
            var id = [];
            id = gijunk.split("|");
            tname = id[1];
        }
        var tooltip = null;

        if (this.m_Config && this.m_Config.TrackConfig) {
            var tcfgs = this.m_Config.TrackConfig;
            Ext.each(tcfgs, function(cfg) {
                //if (cfg.Name == tname) {
                if (cfg.name == tname) {
                    //if (cfg.Help) tooltip = cfg.Help;
                    if (cfg.help) tooltip = cfg.help;
                    return false;
                }
            });
        }

        return tooltip;
    },

    notifyViewLoaded: function(view) {
        if(view != this.m_Panorama) {
            this.fireEvent('graphical_image_loaded', view);
        } else {
            this.fireEvent('panorama_image_loaded', view);
        }
        if (this.m_InitialLoading) {
            var loading_finished = true;
            this.forEachView(function(v) { loading_finished = loading_finished && !v.isLoading(); });
            if (loading_finished) {
                this.m_InitialLoading = false;

                 this.loadSearchPatternData();
                 if (!this.hasMarkers()) {
                    for (var m = 0; m != this.m_Markers.length; m++) {
                        this.addMarkerByData(this.m_Markers[m]);
                    }
                } // create markers from URL
                if (this.m_ViewMarkers) { // show marker info dialog?
                    this.m_ViewMarkers = false;
                    this.showMarkersDlg( view != this.m_Panorama ? view : null );
                }
                if (this.m_ViewSearch && this.m_ViewSearch.length > 0) {
                    var search_term = this.m_ViewSearch;
                    this.forEachView(function(v) {
                        var gotoBox = v.m_View.down('#gotoBox');
                        if (gotoBox) gotoBox.setValue(search_term);
                    });
                    this.gotoAndSearch(this.m_ViewSearch);
                }


            }
            // we need to reload the panorama image in case if the browser added vscroll bar
            if (this.m_Panorama && this.m_Panorama.getWidth() != this.m_Panorama.getScreenWidth()) {
                this.m_Panorama.loadImage();
            }
        }
    },

    moveTo: function(from, len){
        this.forEachView(function(v) {
            v.moveTo(from, len, {from_ui: false});
        });
    },

    addView: function(cfg) {
        var view  = this.m_Panel.add(cfg);
        this.m_Panel.updateLayout();
        return view;
    },

    updateLayout: function() {
        if (this.m_Panel)
            this.m_Panel.updateLayout();
    },

    viewIsClosing: function(view) {
    },

    countActiveSeqPanels: function() {
        var count = 0;
        for (var i = 0; i < this.m_Views.length; i++) {
            if (this.m_Views[i] && this.m_Views[i].m_Type == 'graphical')
                count++;
        }
        return count;
    },

    removeView: function(view) {
        for (var i = 0; i < this.m_Views.length; i++) {
            var v = this.m_Views[i];
            if (v && v.m_Idx == view.m_Idx) {
                this.m_Views[i] = null;
            }
        }
        this.reCreateReflections();
    },

    createPanorama: function() {
        this.m_Panorama = new SeqView.Panorama(this);
        this.m_Views.push(this.m_Panorama);
        this.m_Panorama.loadImage();
    },

    loadPanoramaImage: function() {
        if (this.m_Panel && this.m_Panorama) {
            this.m_Panorama.loadImage();
        }
    },

    registerView: function(view) {
        this.m_Views.push(view);
        view.createPanel();
        this.reCreateReflections();
    },

    createTextView: function(view) {
        // when the TextView dialog (view) closes it will set this.m_TextView to null
        if (!this.m_TextView)
            this.m_TextView = new SeqView.TextView(this);
        this.m_TextView.openDlg(view);
    },

    getPanoramaHeight: function() {
        return this.m_Panorama? this.m_Panorama.getHeight() : 0;
    },

    getPanoramaWidth: function() {
        return this.m_Panorama? this.m_Panorama.getWidth() : 0;
    },

    updateLocator: function(view) {
        if (view.isPanorama() || view.isLoading()) {
            return;
        }

        var locator = view.m_Locator;
        if (locator && this.m_Panorama) {
            if (!this.m_Panorama.isLoading()) {
                var range = view.toSeq();
                locator.setLeft(this.m_Panorama.toPix(range[0]));
                locator.setWidth(this.m_Panorama.toPix(range[2]));
                locator.setHeight(this.getPanoramaHeight());
            } else {
                Ext.defer(this.updateLocator, 50, this, [view]);
            }
        }
    },

    updateViewTitles: function() {
        this.forEachView(function(view) { view.updateTitle(); });
    },

    findView: function(idx) {
        for (var i = 0; i < this.m_Views.length; i++) {
            var view = this.m_Views[i];
            if (view && view.m_Idx == idx) {
                return view;
            }
        }
        return null;
    },

    forEachView: function(fn, scope) {
        var array = this.m_Views;
        for(var i = 0, len = array.length; i < len; i++){
            var v = array[i];
            if (v) {
                if (fn.call(scope || v, v, i, array) === false) { return (i+1); }
            }
        }
    },

    getFlip: function() {
        return this.m_Flip;
    },

    getOrigin: function() {
        return this.m_Origin;
    },

    setFlip: function(flip, view) {
        this.m_Flip = flip;
        var view_idx = -1;
        if (typeof view != "undefined") {
            view_idx = view.m_Idx;
        }
        this.forEachView(function(view) {
            if (view.isGraphic()  &&  view.m_Idx != view_idx) {
                view.setFlipLocal(flip);
            }
        });
        this.fireEvent('strand_changed', this, view);
    },

    decodeRange: function(range) {
        var from = NCBIGBUtils.stringToNum(range[0])-1;
        var to = 0;
        if (!from) from = 0;
        if (from >= this.m_SeqLength) from = 0;
        if (range[1] == 'zs') {
            var len = Math.floor((this.m_Panel.body.getWidth() - 10) * SeqView.MinBpp);
            to = from+len-1;
        } else {
            to = NCBIGBUtils.stringToNum(range[1])-1;
        }
        if (!to) to = this.m_SeqLength;
        to = Math.min(to, this.m_SeqLength - 1);
        return [from, to];
    },

    // Keep pos or range regex in one place
    splitPosOrRange: function(s) {
        s = s.replace(/[, ]/g, '');
        return s.match(/^([-+]?\d+(?:\.\d+)?[km]?)(?:(-|to|\.\.+|:|\/|_)([-+]?\d+(?:\.\d+)?[km]?))?$/i);
    },

    // s - position or range to handle
    // options - Object with parameters
    //     allow_equal - should we allow both ends of range to be equal, false
    //     ask_user    - if position is ambiguous, can we ask user, false
    //     success     - callback for successfully parsed position
    //         pos_range - array with 1 or two elements representing pos or range
    //                     in sequence coordinates
    //         options   - options passed to handlePos
    //     failure       - callback for reporting errors
    //         error message
    //         options   - options passed to handlePos
    //     scope       - object to use as 'this' for callbacks
    handlePos: function(s, options) {
        var parts = this.splitPosOrRange(s);
        if (parts == null) {
            var hgvs = this.isHGVSExpression(s);
            var snip = s.match(/^([rs]s)([0-9]{3,})(?::.+)?$/);
            var vardb = s.match(/^([en]s)(td|v|sv)([0-9]+)(?::.+)?$/);
            if (hgvs) {
                this.handleHGVS(s, options);
            } else if (snip) {
                this.handleVariations(snip[1]+snip[2], options);
            } else if (vardb) {
                this.handleVariations(vardb[1]+vardb[2] + vardb[3], options);
            } else if (options.failure) {
                options.failure.call(options.scope, "Invalid position/range", options);
            }
            return;
        }
        var pos0 = parts[1];
        var sep = parts[2];
        var pos1 = parts[3];
        if (this.m_Origin && !this.isExplicitPosOrRange(pos0, sep, pos1)) {
            if (options.ask_user) {
                // Construct choices
                var is_range = pos1 && pos1.length > 0;
                var ambiguousItems = [];
                var lim0 = 2, lim1 = 2;
                var num0 = NCBIGBUtils.stringToNum(pos0);
                var num1 = NCBIGBUtils.stringToNum(pos1);
                var fc0 = pos0.charAt(0);
                if (fc0 == '+' || fc0 == '-') lim0 = 1;
                var fc1;
                if (is_range && sep !== '/') fc1 = pos1.charAt(0);
                if (fc1 == null || (fc1 == '+' || fc1 == '-')) lim1 = 1;
                var prefix = SeqView.base_url;
                for (var i = 0; i < lim0; ++i) {
                    for (var j = 0; j < lim1; ++j) {
                        var item;
                        var flipped = this.getFlip() ? "_flip" : '';
                        if (is_range) {
                            var n0 = (1-i*2)*num0;
                            var n1 = (1-j*2)*num1;
                            if (n1 < n0) {
                                var t = n0; n0 = n1; n1 = t;
                            }
                            var imageSel = 'mp';
                            if (n0 >= 0 && n1 >= 0) {
                                imageSel = "pp";
                            } else if (n0 < 0 && n1 < 0) {
                                imageSel = "mm";
                            } else if (Math.abs(n0) < Math.abs(n1)) {
                                imageSel = "pm";
                            } 
                            if (n0 > 0) n0 = '+' + n0;
                            if (n1 > 0) n1 = '+' + n1;
                            var fullPos = n0 + ' .. ' + n1;
                            var image = '<img src="' + prefix + 'images/choice_range_' + imageSel + flipped + '.png"/>';
                            item = {boxLabel: image + '&nbsp;' + fullPos, name: 'pos', inputValue: fullPos, checked: (i == 0 && j == 0)};
                        } else {
                            var n0 = (1-i*2)*num0;
                            var imageSel = 'm';
                            if (n0 > 0) {
                                n0 = '+' + n0;
                                imageSel = 'p';
                            }
                            var image = '<img src="' + prefix + 'images/choice_pos_' + imageSel + flipped + '.png"/>';
                            item = {boxLabel: image + '&nbsp;' + n0, name: 'pos', inputValue: n0, checked: (i == 0 && j == 0)};
                        }
                        ambiguousItems.push(item);
                    }
                }
                var title = "Ambiguous " + (is_range ? "range" : "position");
                var subtitle = "Choose the correct " + (is_range ? "range" : "position");
                var dlg = new Ext.Window({
                    app: this,
                    modal: true,
                    title: title,
                    width: 400,
                    layout: 'fit',
                    items: [{
                        xtype: 'form',
                        bodyPadding: '5px',
                        items: [
                            {xtype: 'fieldset', title: subtitle, layout: 'form',
                                items: [{ xtype: 'radiogroup', columns: 1, items: ambiguousItems }]},
                            {xtype: 'displayfield',
                                value: 'To avoid ambiguity, the range coordinate values must use explicit + and - signs'
                                + ' when the origin has been reset.'
                        }]
                    }],
                    buttons: [
                        {text: 'OK', scope: this,
                            handler: function() {
                                var values = dlg.items.items[0].getForm().getValues();
                                this.handlePos.call(this, values['pos'], options);
                                dlg.close();
                            }},
                        {text: 'Cancel', handler: function() { dlg.close(); }}
                    ]
                });
//                if (document.body.clientHeight > dlg.height) dlg.renderTo = this.m_DivId;
                dlg.show();
            } else {
                if (options.failure) {
                    options.failure.call(options.scope, "Explicit position/range required", options);
                }
            }
        } else {
            var beg_pos = this.convertRelativePosition(pos0);
            if (pos1 && pos1.length > 0) {
                // range is specified
                var end_pos;
                if (sep == "/") {
                    var pad = NCBIGBUtils.stringToNum(pos1);
                    end_pos = beg_pos + pad;
                    beg_pos -= pad;
                } else {
                    end_pos = this.convertRelativePosition(pos1);
                    if (beg_pos > end_pos) {
                        var t = beg_pos;
                        beg_pos = end_pos;
                        end_pos = t;
                    }
                }
                if (isNaN(beg_pos) || isNaN(end_pos) || beg_pos < 0 || end_pos < 0 ||
                    (options.allow_equal ? beg_pos > end_pos : beg_pos >= end_pos)  ||
                    beg_pos >= this.m_SeqLength || end_pos >= this.m_SeqLength)
                {
                    if (options.failure) {
                        var from_pos = this.posToLocal(0);
                        var to_pos = this.posToLocal(this.m_SeqLength - 1);
                        if (to_pos < from_pos) {
                            var temp = to_pos;
                            to_pos = from_pos;
                            from_pos = temp;
                        }
                        var msg = 'Invalid range: Sequence positions should be from ' + from_pos + ' to ' + to_pos;
                        options.failure.call(options.scope, msg, options);
                    }
                } else {
                    if (options.success) {
                        options.success.call(options.scope, [beg_pos, end_pos], options);
                    }
                }
            } else {
                // single position
                if (isNaN(beg_pos) || beg_pos < 0 || beg_pos >= this.m_SeqLength) {
                    if (options.failure) {
                        var from_pos = this.posToLocal(0);
                        var to_pos = this.posToLocal(this.m_SeqLength - 1);
                        if (to_pos < from_pos) {
                            var temp = to_pos;
                            to_pos = from_pos;
                            from_pos = temp;
                        }
                        var msg = 'Invalid position: Sequence position should be from ' + from_pos + ' to ' + to_pos;
                        options.failure.call(options.scope, msg, options);
                    }
                } else {
                    if (options.success) {
                        options.success.call(options.scope, [beg_pos], options);
                    }
                }
            }
        }
    },

    // term - HGVS position to parse
    // options - Object with parameters
    //     allow_equal - should we allow both ends of range to be equal, false
    //     success     - callback for successfully parsed position
    //         pos_range - array with 1 or two elements representing pos or range
    //                     in sequence coordinates
    //         options   - options passed to handlePos
    //     failure       - callback for reporting errors
    //         error message
    //         options   - options passed to handlePos
    //     scope       - object to use as 'this' for callbacks
    handleHGVS: function(term, options) {
        var view = options ? options.view : null;
        var params = {term: term, id: this.GI, type: "hgvs"};
        if (this.m_Key && this.m_Key.length > 0)
            params.key = this.m_Key;
        var tracks = this.getActTracks();
        if (tracks)
            params.tracks = tracks;

        return [this.AjaxRequest({
            url: this.m_CGIs.FeatSearch,
            data: params,
            context: this,
            success: function(data, text, res) {
                var from_cgi = SeqView.decode(data);
                if (from_cgi.total) {
                    var feat = from_cgi.features[0];
                    var pos = feat.from;
                    var pos_range = [pos];
                    if (feat.to != pos) {
                        pos_range.push(feat.to);
                    }
                    var label = feat.label;
                    if (options && options.success) {
                        options.success.call(this, pos_range, options, label);
                    } else {
                        this.setPositionalMarker(pos_range, label);
                        if (view) {
                            view.gotoPosRange([pos], true, {from_ui: true} );
                        } else {
                            this.forEachView(function(v) {
                                if (v.isGraphic())
                                    v.gotoPosRange([pos], true, {from_ui: true});
                            });
                        }
                    }
                } else {
                    var msg = from_cgi.errmsg || (term + " not found on the sequence");
                    if (options && options.failure) {
                        options.failure.call(this, msg, options);
                    } else {
                        Ext.MessageBox.alert('Search ' + (from_cgi.errmsg ? 'Error' : 'results'), msg);
                    }
                }
            },
            error: function(data, text, res) {
                var msg = "Server error in parsing/searching request " + term;
                if (options && options.failure) {
                    options.failure.call(this, msg, options);
                }
            }
        })];
    },

    handleVariations: function(term, options) {
        var rqNum = 1;
        var view = options ? options.view : null;
        var app = this;
        var extSearch = false,
            found = 0;
        var pos_range, params;
        var completeTask = function() {
            if (--rqNum > 0) return;
            if (!found && !extSearch) { //searching in not displayed tracks
                rqNum = 0;
                extSearch = true;
                Ext.each(app.m_Config.TrackConfig, function(trk, idx) {
                    if (trk.shown || trk.category.name != 'Variation') return true;
                    params.track = '[id:' + trk.id + ',key:' + trk.key + ',filter:' + trk.filter + ']';
                    params.pos = params.from;
                    rqNum++;
                    app.Request({url: app.m_CGIs.Graphic, data: params,
                        success: function(res) { processResponse(res, idx); },
                        error: completeTask});
                });
                return;                
            }
            if (found && extSearch) app.fireEvent('configuration_changed', app);
            if (view) {
                view.gotoPosRange(pos_range, true, {from_ui: true} );
            } else {
                app.forEachView(function(v) {
                   if (v.isGraphic()) v.gotoPosRange(pos_range, true, {from_ui: true});
                });
            }
        }
        var processResponse = function(res, idx) {
            if (res.error_message) {
                app.unmaskGraphViews(view);
                Ext.MessageBox.alert('Search Error', res.error_message);
                return;
            }
            if (res.new_pos == pos_range[0]) {
                app.m_Config.TrackConfig[idx].shown = true;
                found++;
                if (!extSearch) rqNum = 1;
                //console.log (rqNum + '. ' + app.m_Config.TrackConfig[idx].id + ': ' + res.new_pos + ' ' + app.m_Config.TrackConfig[idx].display_name);
            } else {
                var jst = res.job_status;
                if (jst == 'submitted' || jst == 'running' || jst == 'pending') {
                    Ext.defer(SeqView.App.simpleAjaxRequest, 1000 * (rqNum % 5), app,
                        [{url: app.m_CGIs.Graphic, data: {job_key: res.job_id},
                        success: function(res) { processResponse(res, idx);},
                        error: completeTask}]);
                    return;
                } 
            }
            completeTask();
        }

        var params = {term: term, id: this.GI, type: "snp"};
        if (this.m_Key) params.key = this.m_Key;
        app.maskGraphViews(view, 'Searching...');
        return [this.AjaxRequest({
            url: this.m_CGIs.FeatSearch,
            data: params,
            context: this,
            success: function(data, text, res) {
                var from_cgi = SeqView.decode(data);
                if (from_cgi.total) {
                    var snp_found = false;
                    for (var i = 0; i < from_cgi.features.length; i++) {
                        var feat = from_cgi.features[i];
                        if (feat.seqId != this.GI)
                            continue;
                        snp_found = true;
                        var pos = feat.from;
                        pos_range = [pos];
                        if (feat.to != pos) {
                            pos_range.push(feat.to);
                        }
                        var label = feat.label;
                        if (options && options.success) {
                            options.success.call(this, pos_range, options, label);
                        } else {
                            this.setPositionalMarker(pos_range, label);
                            params = {id: app.GI, navi: 1, len: 10, width: 100, from: pos + (!pos ? 1 : -1), dir: (!pos) ? 'prev' : 'next'};
                            // searching in displayed tracks
                            Ext.each(app.m_Config.TrackConfig, function(trk, idx) {
                                if (!trk.shown || trk.category.name != 'Variation') return true;
                                params.track = '[id:' + trk.id + ',key:' + trk.key + ',filter:' + trk.filter + ']';
                                params.pos = params.from;
                                rqNum++;
                                app.AjaxRequest({url: app.m_CGIs.Graphic, data: params,
                                    success: function(res) { processResponse(res, idx);},
                                    error: completeTask});
                            });
                            completeTask();
                        }
                        break;
                    }
                    // it's a valid snp but it's not annotated on this sequence
                    // so we suggest exploring dbSNP page
                    if (!snp_found) {
                        app.unmaskGraphViews(view);
                        var snip = term.match(/^([rs]s)([0-9]{3,})$/);
                        var snip_url = SeqView.webNCBI + 'snp/' + snip[2];
                        Ext.Msg.alert("SNP search", 'Valid SNP is not found on this sequence, please see <a target="_blank" href="'+snip_url+'">dbSNP page for '+term+'</a>');
                    }
                } else {
                    var msg = term + " not found on the sequence";
                    if (options && options.failure) {
                        options.failure.call(this, msg, options);
                    }
                }
            },
            error: function(data, text, res) {
                var msg = "Server error in parsing/searching request " + term;
                if (options && options.failure) {
                    options.failure.call(this, msg, options);
                }
            }
        })];
    },

// parsePosOrRange:
// returns an array of one or two elements, depending on input, which
//    represent a valid position or range on sequence. If there is no valid
//    position/range can be parsed, returns undefined. All input is in local
//    coordinate system - 1-based with flip and origin taken into account.
//    Output is 0-based, in sequence coordinates.

    parsePosOrRange: function(s, allow_equal) {
        s = s.replace(/[, ]/g, '');
        var parts = this.splitPosOrRange(s);
        if (parts == null) {
            return;
        }
        var part0 = parts[1];
        var sep = parts[2];
        var part1 = parts[3];
        var beg_pos = this.convertRelativePosition(part0);
        if (part1 && part1.length > 0) {
            // range is specified
            var end_pos;
            if (sep == "/") {
                var pad = NCBIGBUtils.stringToNum(part1);
                end_pos = beg_pos + pad;
                beg_pos -= pad;
            } else {
                end_pos = this.convertRelativePosition(part1);
                if (this.getFlip()) {
                    var t = beg_pos;
                    beg_pos = end_pos;
                    end_pos = t;
                }
            }
            if (isNaN(beg_pos) || isNaN(end_pos) || beg_pos < 0 || end_pos < 0 ||
                (allow_equal ? beg_pos > end_pos : beg_pos >= end_pos)  ||
                beg_pos >= this.m_SeqLength || end_pos >= this.m_SeqLength)
            {
                return;
            }
            return [beg_pos, end_pos]
        } else {
            // single position
            if (isNaN(beg_pos) || beg_pos < 0 || beg_pos >= this.m_SeqLength)
                return;
            return [beg_pos];
        }
    },


    isPosOrRange: function(s) {
        return this.splitPosOrRange(s) != null;
    },

//////////////////////////////////////////////////////////////////////////
// isExplicitPosOrRange:
//     pos0, sep, pos1 - parts of a range split by pattern
//     returns true if numerals have a sign. e.g. -12:+24, or +35..-20
    isExplicitPosOrRange: function(pos0, sep, pos1) {
        var fc0 = pos0.charAt(0);
        var fc1;
        if (pos1 && pos1.length > 0 && sep !== '/') fc1 = pos1.charAt(0);
        return !((fc0 != '+' && fc0 != '-') || (fc1 && fc1 != '+' && fc1 != '-'));
    },

    convertRelativePosition: function(pos_str) {
        var pos = NCBIGBUtils.stringToNum(pos_str);
        if (!isNaN(pos))
            return this.posToGlobal(pos);
        return pos;
    },

//////////////////////////////////////////////////////////////////////////
// posToLocal: convert backend 0-based global coordinate to
//     local, taking into account strand and origin
    posToLocal: function(pos, flip) {
        pos -= this.m_Origin;
        if (pos >= 0) pos += 1;
        return pos;
    },

//////////////////////////////////////////////////////////////////////////
// posToLocalDisplay: convert backend 0-based global coordinate to
//     local, taking into account origin, returns number with sign if
//     origin is set, so for positive numbers it will return a string
    posToLocalDisplay: function(pos, flip) {
        pos -= this.m_Origin;
        if (pos >= 0) {
            pos += 1;
            if (this.m_Origin) pos = '+' + pos;
        }
        return pos;
    },

//////////////////////////////////////////////////////////////////////////
// posToGlobal: convert 1-based coordinate relative to origin and strand
//     to global 0-based backend coordinate
    posToGlobal: function(pos) {
        if (pos > 0) pos -= 1;
        return pos + this.m_Origin;
    },

    newMarkerDlg: function(view, x_pos) {
        if (!this.m_MarkersInfo) {
            this.m_MarkersInfo = new SeqView.MarkersInfo(this);
        }
        this.m_MarkersInfo.MarkerDlg(view || this, x_pos);
    },

    showMarkersDlg: function( view ){
        if( !this.m_MarkersInfo ){
            this.m_MarkersInfo = new SeqView.MarkersInfo( this );
        }
        this.m_MarkersInfo.showDlg( view );
    },

    updateMarkersSize: function(view) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.updateMarkersSize(view);
        }
    },

    hasMarkers: function() {
        return this.m_MarkersInfo && this.m_MarkersInfo.hasMarkers();
    },

    getMarkersInfo: function() {
        return this.m_MarkersInfo;
    },

    updateMarkersPos: function(view) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.updateMarkersPos(view);
        }
    },

    forEachMarker: function(fn, scope) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.forEachMarker(fn,scope);
        }
    },

    findMarker: function(m_name) {
        if (this.m_MarkersInfo) {
            return this.m_MarkersInfo.findMarker(m_name);
        }
        return null;
    },

    addMarker: function(marker_data) {
        if (!this.m_MarkersInfo) {
            this.m_MarkersInfo = new SeqView.MarkersInfo(this);

        }
        this.m_MarkersInfo.addMarker(marker_data);
    },

    addMarkerByData: function(marker_data) {
        if (!this.m_MarkersInfo) {
            this.m_MarkersInfo = new SeqView.MarkersInfo(this);
        }
        this.m_MarkersInfo.addMarkerByData(marker_data);
    },

    scrollMarkers: function(view, delta) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.scrollMarkers(view,delta);
        }
    },


    updateReflections: function() {
        if( this.m_Reflections ){
            this.m_Reflections.updateAll();
        }
   },

    scrollReflections: function(view, delta) {
        if (this.m_Reflections) {
            this.m_Reflections.scrollPix(view,delta);
        }
        //this.reCreateReflections();
   },

    reCreateReflections: function() {
        return; // temporary disabled
        if (!this.m_Reflections)
            this.m_Reflections = new SeqView.ReflectionCont(this);
        this.m_Reflections.reCreate();
    },

    reloadAllViews: function(options) {
        this.forEachView(function(view) { view.refresh(options); });
    },

    showOriginDlg: function(pos) {
        Ext.MessageBox.prompt('Set Origin', 'Please enter new origin:', function(btn, text) {
            if (btn!='ok'  || text.length==0) return;
            var position = 0;
            var bad_pos = !NCBIGBUtils.isNumeric(text);
            if (!bad_pos) {
                position = NCBIGBUtils.stringToNum(text);
                if (position > this.m_SeqLength) {
                    bad_pos = true;
                }
            }
            if (bad_pos) {
                Ext.MessageBox.show({title: 'Set Origin',msg: 'Invalid sequence origin.',
                                     buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                return;
            }
            position -= 1;
            if (this.m_Origin == position)
                return;
            this.m_Origin = position;
            if (this.m_Origin < 0) this.m_Origin = 0;
            // reload all views
            this.reloadAllViews();
            this.fireEvent('origin_changed', this);
        }, this, false, pos ? pos+1 : this.m_Origin+1);
    },

    setOrigin: function(view, x_pos) {
        var pix = x_pos - view.m_ScrollPix;
        var seq_pos = view.pix2Seq(pix);
        this.showOriginDlg(seq_pos);
    },

    clearOrigin: function() {
        Ext.MessageBox.confirm('Confirm', 'Reset Sequence Origin?', function(btn) {
            if (btn!='yes' || this.m_Origin == 0) return;
            this.m_Origin = 0;
            // reload all views
            this.reloadAllViews();
            this.fireEvent('origin_changed', this);
        }, this);
    },

    loadExtraData: function() {
        this.load_params = this.load_params || {};
        if (!this.m_userData.length) {
            if (this.GI) this.createAndLoadMainPanel();
            else this.infoMissed();
            return;
        }

        if (typeof this.m_AssmContext == 'undefined') {
            this.m_AssmContext = '';
            if (this.GI && this.GI.indexOf("lcl|") !=0) {
                var options = { context: this,
                    url: this.m_CGIs.Config + '?notrackconfig=true&id=' + this.GI,
                    success: function(data, txt, rq) {
                        if (data.job_status) {
                            var st = data.job_status;
                            if (st == 'submitted' || st == 'running' || st == 'pending') {
                                options.url = this.m_CGIs.Config + '?job_key=' + data.job_id;
                                Ext.defer(SeqView.App.simpleAjaxRequest, 1000, this, [options]);
                                return;
                            }
                       }
                       if (data.SeqInfo) this.m_AssmContext = data.SeqInfo.assm_context;
                       this.loadExtraData();
                    },
                    error: this.loadExtraData
                };
                this.AjaxRequest(options);
                return;
            }
        }

        var config = {assm_acc: this.m_AssmContext};//, data_action:'uploading', ret_prj_key: 'true' };
        if (this.GI) config['accession'] = this.GI;
        if (!(config.prj_key = this.m_Key || this.load_params.key)) delete config.prj_key;
        var data_item = this.m_userData.shift();
        if (data_item[0] == 'rkey') {// Fetch long RID with filters from NetCache key
            var params = {data_action: 'downloading', format: 'rids', key: data_item[1], fmt: "text/plain"};
            var processResponse = function(data) {
                    if (data.statusText == 'OK') this.m_userData.unshift(['rid', data.responseText]);
                    else  this.m_userTrackNames.shift(); // Keep name array in sync
                    this.loadExtraData();
            }
            this.AjaxRequest({url: this.m_CGIs.NetCache, data: params, context: this,
                success: processResponse,
                error: processResponse});
        } else {
            config.track_name = this.m_userTrackNames.shift();
            if (config.track_name) delete config.track_name;
            switch (data_item[0]) {
                case 'rid': config.blast = {rid: data_item[1], link_related_hits: this.m_FindCompartments};
                    this.load_params.rid_loaded = true; // Mark that we loaded RID
                    break;
                case 'data':
                    config.data = data_item[1];
                    config.file_format = 'asn text';
                    break;
                case 'url_reload': config.check_cs = true; // no break;
                case 'url': config.dataURL = data_item[1]; break;
            }
            var msg = 'Uploading your data';
            if (this.m_AssmContext) {
                msg += " on assembly " + this.m_AssmContext;
            }
            msg += ', please wait...'; 
            this.showMessage(msg);
            this.m_Panel.items.items[0].getEl().mask('Uploading \"' + data_item[0] + '\"');
            var app = this;
            app.consError = console.error;
            console.error = function() {}
            try {
                var uploaderUUD = new UUD.FileUploader(config);
                var promise = uploaderUUD.getPromise();
                promise.fail(function() {
                    console.error = app.consError || console.error; delete app.consError;
                    app.m_Panel.items.items[0].getEl().unmask();
                    var errMsg = this.getErrors();
                    app.showMessage(errMsg, true);
                    var fbText = 'SViewer initial parameters: ' + app.m_AllViewParams + '\nError message: ' + this.getErrors();
                    Ext.Msg.show({title: 'User data loading error',
                        msg: errMsg,
                        buttons: Ext.Msg.YESNOCANCEL,
                        buttonText: {yes: 'Continue', no: 'Feedback'},
                        icon: Ext.MessageBox.WARNING,
                        fn: function(btn) {
                            switch (btn) {
                                case 'no': app.showFeedbackDlg('2-Feedback', fbText, document.location.href); // "break" is missed to continue data loading
                                case 'yes': app.loadExtraData(); break;
                            }
                        }});
                });
                promise.done(function(tlist, dkey) {
                    console.error = app.consError || console.error; delete app.consError;
                    if (app.load_params.rid_loaded && !app.GI) app.initRID(this.getMetadata());
                    app.addKey(dkey);
                    var tracks = this.getTracks();
                    if (tracks) Ext.each(tracks, function() { app.addUploadedTrackID(this); });
                    var cookie_name = SeqView.Cookies.UserDataCookieName;
                    if (cookie_name === SeqView.Cookies.UserDataCookieNameBase) {
                        cookie_name += (app.m_AppName ? ('-' + app.m_AppName) : '');
                    }
                    SeqView.SessionData.set(cookie_name, app.m_Key);//save the key into a session cookie
                    var msg = 'Data uploaded';
                    if (app.m_AssmContext) 
                        msg += " on assembly " + app.m_AssmContext;
                    app.showMessage(msg);
                    app.loadExtraData();
                });
                var currTask = 'Uploading your data'
                promise.progress(function(progress) {
                    var task = progress.current_task;
                    if (task == "" || task == 'pending' || task == currTask) return;
                    app.showMessage(currTask = task);
                });
                uploaderUUD.upload();
            } catch(e) { app.showMessage('Unable to upload data: ' + e.message, true); } 
        }
    },
    
    addUploadedTrackID: function(track) {
        var id = track.id || track.GetTMSId();
        this.uploadedIDs = this.uploadedIDs || [];
        var obj = this.uploadedIDs.filter(function ( obj ) {
            return obj.id === id;
        })[0];
        if (obj === undefined)
            this.uploadedIDs.push(track);
        //if (this.uploadedIDs.indexOf(id) == -1) this.uploadedIDs.push(id);
    },

    showMessage: function(msg, errorFlag) {
        var msg_field = Ext.getCmp((errorFlag ? 'sv-uplerr' : 'sv-uplmsg') + this.m_Idx);
        if (msg_field) {
            msg_field.show();
            msg_field.update(msg);
        }
    },

    maskGraphViews: function(view, msg, msgCls) {
        var mask = function(v) { if (v.isGraphic()) Ext.get(v.m_DivId).mask(msg, msgCls); }
        if (view) mask(view); else this.forEachView(mask);
    },
    unmaskGraphViews: function(view) {
        var unmask = function(v) { if (v.isGraphic()) Ext.get(v.m_DivId).unmask(); }
        if (view) unmask(view); else this.forEachView(unmask);
    },

    addKey: function(key) {
        if (typeof key == 'undefined') return;
        if (!this.m_Key || this.m_Key.length == 0) this.m_Key = key;
        else
            if (this.m_Key.search(key) == -1) this.m_Key += ',' + key;
    },

    initRID: function(from_cgi) {
        if (from_cgi['blast_query'] === undefined) return; // Can't do anything
        var blast_query_list = from_cgi['blast_query'];
        var blast_query = {};
        for (var i = 0; i < blast_query_list.length; i++) {
            var kv = blast_query_list[i];
            blast_query[kv.key] = kv.value;
        }
        this.GI = blast_query['id'];
        if (!this.m_QueryRange) {
            var beg = blast_query['beg'];
            var end = blast_query['end'];
            var querybeg = blast_query['querybeg'];
            var queryend = blast_query['queryend'];
            if (querybeg !== undefined  &&  queryend !== undefined) {
                // Convert to numbers
                querybeg = querybeg - 0;
                queryend = queryend - 0;

                // Adjust viewrange, add special marker to show query if requested
                beg = beg !== undefined ? Math.min(beg, querybeg) : querybeg;
                end = end !== undefined ? Math.min(end, queryend) : queryend;
                if (this.m_ShowQuery) {
                     this.m_Markers.push([ [querybeg, queryend], "Query", SeqView.MarkerFlags.SystemLock | SeqView.MarkerFlags.Hollow, "red"]);
                }
            }
            if (beg !== undefined  &&  end !== undefined) {
                var overhang = 0.15 * (end - beg);
                beg = Math.max(0, beg - overhang);
                end += overhang;
                beg = Math.round(beg);
                end = Math.round(end);
                this.m_ViewRanges = [ beg + ':' + end ];
            }
        }
    },

    onSyncAlignView: function(view) {
        for (var i = 0; i < this.m_Views.length; ++i) {
            aview = this.m_Views[i];
            if (aview && aview.isAlignment() ) {
                aview.loadImage(view.m_VisFromSeq, view.m_VisLenSeq);
                break;
            }
        }
    },

    getGraphicViews: function() {
        var views = []
        this.forEachView(function(view) {
            if (view.isGraphic())
                views.push(view);
        });
        return views;
    },

    AjaxRequest: function(cfg){
/*
        if( !this.m_AjaxTrans ) this.m_AjaxTrans = [];

        cfg.callback = cfg.callback || function(arg) {};
        cfg.params = cfg.params || {};
*/
        Ext.applyIf (cfg.data, {appname: this.m_AppName});

/*
        var trans = this.m_AjaxTrans;

        var prefunction = function(opts, succsess, res) {
            if( opts.params.transId ){
                trans.remove( opts.params.transId );
                delete opts.params.transId;
            }
        }
        cfg.callback = Ext.Function.createInterceptor(cfg.callback, prefunction);
*/
        var transId = SeqView.App.simpleAjaxRequest(cfg);
/*        if( transId ){
            cfg.params.transId = transId;
            trans.push( transId );
            return transId;
        }*/
    },

    getCustomToolTipTools: function(selection) {
        if (this.m_CustomSelectionHandler)
            return this.m_CustomSelectionHandler.getToolTipTools(selection);
        return null;
    },
    setCustomSelectionHandler: function(handler) {
        this.m_CustomSelectionHandler = handler;
    },

    addCustomFeatureFlags: function(cfg) {
        if (this.m_CustomSelectionHandler) {
            this.m_CustomSelectionHandler.addCustomFeatureFlags(cfg);
        }
    },

    setTooltipPreprocessor: function(callback) {
        if (typeof callback === 'function') this.m_preprocessorTT = callback;
        else delete this.m_preprocessorTT;
    },

    watchdogStart: function(url, job_id, params) {
        this.watchdogStop();
        this.m_WatchUrl = url;
        this.m_WatchJobId = job_id;
        this.m_WatchJobParams = params;
        this.m_WatchTimeoutId = Ext.defer(this.watchdogReport, 90000, this, ['No response']);
    },

    watchdogStop: function() {
        if (this.m_WatchTimeoutId ){
            clearTimeout(this.m_WatchTimeoutId);
            this.m_WatchTimeoutId = null;
        }
    },

    watchdogReport: function(reason) {
        this.watchdogStop();
        this.AjaxRequest({
            url: this.m_CGIs.Watchdog,
            data: {
                reason: reason,
                requrl: this.m_WatchUrl,
                reqparams: this.m_WatchJobParams,
                jobid: this.m_WatchJobId
            }
        });
    }

});

Ext.apply(SeqView.App, {

    findAppByDivId: function(div_id) {
        for (var i = 0; i < SeqView.m_Apps.length; ++i) {
            if (SeqView.m_Apps[i].m_DivId == div_id) return SeqView.m_Apps[i];
        }
        return null;
    },

    findAppByIndex: function(app_idx) {
        for (var i = 0; i < SeqView.m_Apps.length; ++i) {
            if (SeqView.m_Apps[i].m_Idx == app_idx) return SeqView.m_Apps[i];
        }
        return null;
    },

    getApps: function() { return SeqView.m_Apps; },

    simpleAjaxRequest: function(cfg) {
        cfg.xhrFields = {withCredentials: true};
        Ext.applyIf(cfg, {type: 'POST'});
        cfg.crossDomain = true;
        Ext.applyIf(cfg, {dataType: SeqView.jsonType});
        jQuery.support.cors = true;
        try {
            return jQuery.ajax(cfg);
        } catch (e) { return e;};
    },

    showLinkURLDlg: function(app_idx, logarea) {
        var app = this.findAppByIndex(app_idx);
        if (app) { app.showLinkURLDlg(logarea); }
    },

    showFeedbackDlg: function(app_idx, logarea) {
        var app = this.findAppByIndex(app_idx);
        if (app) app.showFeedbackDlg(logarea);
    }
});

// Portal integration function. It is called then Portal needs to know current
// sequence viewer position.
SeqView.PortalSeqGraphicsInfo = function() {
    var params = "from=&to=&itemid=&strand=";
    if(SeqView.App) {
        var apps = SeqView.App.getApps();
        if (apps && apps.length == 1) {
            var app = apps[0];
            var views = app.getGraphicViews();
            if (views && views.length == 1) {
                var view = views[0];
                if (view.m_VisFromSeq >= 0 && view.m_VisLenSeq && view.m_VisLenSeq > 0) {
                    params = 'from='+(view.m_VisFromSeq+1)+'&to='+(view.m_VisFromSeq+view.m_VisLenSeq);
                } else {
                    params = 'from=&to=';
                }
                params += '&itemid=';
                if (app.m_ItemID) {
                    params += app.m_ItemID;
                }
                if (view.canFlip()) {
                    if (view.getFlip()) {
                        params += '&strand=true';
                    } else {
                        params += '&strand=false';
                    }
                } else {
                    params += '&strand=';
                }
            }
        }
    }
    return params;
};

/*  $Id: sviewapp_more.js 38325 2017-04-25 22:08:59Z borodine $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 *
 * File Description:
 *
 */

Ext.apply(SeqView.App.prototype, {

    showLinkURLDlg: function (logarea) {
        if (logarea) this.pingClick(logarea);
        this.getLinkToThisPageURL(null, null, function (linkUrl, longUrl) {
            this.resizeIFrame(400);
            longUrl = longUrl || linkUrl;
            var longParams = longUrl.substr(longUrl.indexOf('?') + 1) + '&appname=no_appname';
            var template = '<iframe id="sv@iframe" width="' + this.m_Views[0].m_Width + '" src="' + SeqView.base_url
                      + 'embedded_iframe.html?iframe=sv@iframe&' + longParams + '" onload="'
                      + 'if(!window._SViFrame){_SViFrame=true;window.addEventListener(\'message\','
                      + 'function(e){if(e.origin==\'https://' + document.domain + '\' && !isNaN(e.data.h))'
                      + 'document.getElementById(e.data.f).height=parseInt(e.data.h);});}">\n</iframe>';
            var embed_div_template = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n'
                      + '<html>\n<head>\n<title>Embedded Sequence Viewer with parameters</title>\n'
                      + '<script type="text/javascript" src="https://www.ncbi.nlm.nih.gov/projects/sviewer/js/sviewer.js"></script>'
                      + '</head>\n<body>'
                      + '<div id="@sviewer_id@" class="SeqViewerApp" data-autoload>\n<a href="?embedded=true&' + longParams + '"></a>\n</div>\n</body>\n</html>';
            var embed_div = "";
            var embed_iframe = "";


            var showEmbedCode = function (uname) {
                if (!uname || uname.length < 1) uname = 'sv_iframe';
                var comp = Ext.getCmp('SV_embed_code');
                embed_iframe = template.replace(/sv@iframe/g, uname);
                embed_div = embed_div_template.replace(/@sviewer_id@/g, uname);
                comp.setValue(embed_iframe);
            }
            var app = this;
            var linkDlg = new Ext.Window({
                app: this,
                layout: 'vbox',
                modal: true,
                title: 'Link To This View',
                width: 700, autoHeight: true,
                constrain: true,
                cls: 'SeqViewerApp',
                listeners: {
                    'beforeshow': function (qt) {
                        NCBIGBUtils.makeTinyURL(linkUrl, function (res) {
                            var uname = '';
                            if (res.id) {
                                Ext.getCmp('the_tinyURL_link').setValue(res.id);
                                var parts = res.id.split('/');
                                uname = parts[parts.length - 1].trim();
                            }
                            showEmbedCode(uname.toLowerCase());
                        }
                        );
                    },
                    'show': function () {
                        if (this.m_Key && this.m_Key.length > 0) {
                            Ext.Msg.show({
                                title: 'Warning',
                                msg: "This page contains a link to the user loaded data " +
                                     "which is kept <br/>inside a temporary storage and will not be available in approximately 1-2 months.",
                                icon: Ext.Msg.WARNING,
                                buttons: Ext.Msg.OK
                            });
                        }
                        if (longUrl && longUrl != linkUrl) {
                            Ext.Msg.show({
                                title: 'Warning',
                                msg: "The tracks configuration string for this link is too long and it has been saved<br/>" +
                                     "into a temporary storage. It will be kept there for 3 months.",
                                icon: Ext.Msg.WARNING,
                                buttons: Ext.Msg.OK
                            });
                        }

                    },
                    scope: this
                },
                resizable: true,
//                resizeHandles: 'ew',
                closeAction: 'destroy',
                plain: true,
                bodyStyle: 'padding:5px;',
//                labelWidth: 140,
                autoHeight: true,
                frame: true,
                labelAlign: 'right',
                items: [
                    { xtype: 'textarea', height: 115, allowBlank: true, width: '100%', flex:1, value: linkUrl, name: 'the_link' },
                    { xtype: 'textfield', id: 'the_tinyURL_link', allowBlank: true, width: '100%', emptyText: 'Loading...', value: '', name: 'the_tinyURL' },
                    { xtype: 'fieldset',  width: "100%", flex:2,
                      autoHeight: true,
                      title: 'Embedding code',
//                      bodyStyle: 'padding:5px',
                      layout: 'vbox', // 'form',
                      items: [
                        { xtype: 'radiogroup', cls: 'x-check-group-alt', width: '100%',
                            items: [
                              { boxLabel: "Embed&nbsp;code&nbsp;for&nbsp;IFRAME", name:'rb-embed', checked: true, inputValue: 'iframe',
                                  listeners: {
                                      change: function(fld, newval, oldval) {
                                          if (!newval) return;
                                          var comp = Ext.getCmp('SV_embed_code');
                                          comp.setValue(embed_iframe);
                                      }}},
                              { boxLabel: "Whole&nbsp;page&nbsp;example", name:'rb-embed', inputValue: 'div',
                                  listeners: {
                                      change: function(fld, newval, oldval) {
                                          if (!newval) return;
                                          var comp = Ext.getCmp('SV_embed_code');
                                          comp.setValue(embed_div);
                                      }}}
                            ]},
                        { xtype: 'textarea', id: 'SV_embed_code', allowBlank: true, flex:1, width: '100%', emptyText: 'Loading...', value: '', name: 'embed_code' }
                    ]}
                ],
                buttons: [{ text: 'Close', handler: function () { linkDlg.close(); app.resizeIFrame(); } }]
            });
            linkDlg.show();
        });
    },


    getLinkToThisPageURL: function (mode, tracks, callback) {
        callback = callback || function (link) { window.open(link); };
        tracks = tracks || SeqView.TM.tracksArrayToString(this.getActiveTracks(true), true, true);
        var appParams = this.getAppParams();
        var link = SeqView.base_url + '?id=' + this.GI;
        var info = this.m_Config.SeqInfo;
        if ((this.m_Portal || mode == 'portal') && mode != 'full' && info.gi) {
            var db = info.acc_type == 'DNA' ? 'nuccore' : 'protein';
            link = SeqView.webNCBI + db + '/' + info.gi + '?report=graph';
        }
        var longUrl = link + '&tracks=' + tracks + appParams;
        if (longUrl.length <= 2016)
           callback.call(this, longUrl); 
        else {
            var params = { data_action: 'uploading', input_form: 'tracks', tracks: tracks };
            this.AjaxRequest({ url: UUD.GetUUDCgiUrl('sv-dl'), data: params, context: this,
                success: function (data, text, res) {
                    callback.call(this, data.success ? link + '&tkey=' + data.key + appParams : longUrl, longUrl);
                },
                error: function (data) { callback.call(this, longUrl); }
            });
        }
    },
    
    
    getAppParams: function() {
        var params = '';
        if (this.m_Origin != 0) params += '&origin=' + this.m_Origin;
        if (this.m_NAA) params += '&naa=' + this.m_NAA;
        if (this.m_SRZ) params += '&srz=' + this.m_SRZ;
        if (this.m_BamPath) params += '&bam_path=' + this.m_BamPath;
        if (this.m_DepthLimit) params += '&depthlimit=' + this.m_DepthLimit;
        if (this.m_SnpFilter) params += '&snp_filter=' + this.m_SnpFilter;
        if (this.m_Key) params += '&key=' + this.m_Key;
        if (this.m_MarkersInfo) params += this.m_MarkersInfo.getMarkersURL();
        if (this.m_ViewerContext) params += '&viewer_context=' + this.m_ViewerContext;
        if (this.m_AssmContext) params += '&assm_context=' + this.m_AssmContext;
        params += (this.m_noGuessAssm ? '&noguess_assm=1' : '');
        // add views
        var view = '&v';
        var color = '&c';
        var slimMode = '&slim';
        var select = '&select';
        this.forEachView(function (v) {
            if (!v.isGraphic()) return;
            if (v.m_UrlFrom && v.m_UrlTo) { // if we have presise positions from URL, use them. Otherwise - calculate from visible range
                view += ',' + v.m_UrlFrom + ':' + v.m_UrlTo;
            } else {
                view += ',' + (v.m_VisFromSeq + 1) + ':' + (v.m_VisFromSeq + v.m_VisLenSeq);
            }
            color += ',' + v.m_Color;
            slimMode += ',' + (0 + v.m_slimMode);
            select += ',' + (v.m_SelectedSig ? v.m_SelectedSig : 'null');
        });

        var options = this.m_Config.Options;
        if (options) {
            Ext.each(['color', 'label', 'decor', 'spacing'], function (attr) {
                var idx = options.controls[attr].indexOf(options['curr_' + attr]);
                if (idx > 0) params += '&' + attr + '=' + idx;
            });
        }
        params += (view.replace(/,/, '=') + color.replace(/,/, '=') + (this.m_Flip ? '&gflip=1' : '')
                 + select.replace(/,/, '=') + slimMode.replace(/,/, '='));

        if (this.m_TextView) params += this.m_TextView.getURL();

        if (this.m_SearchDlg && this.m_SearchDlg.isVisible()) {
            var view = null;
            this.forEachView(function (v) { if (!view && v.isGraphic()) view = v; });
            if (view) {
                var combo = view.m_View.down('#gotoBox');
                var searchTerm = null;
                if (combo) searchTerm = combo.getValue();
                if (searchTerm) params += '&search=' + searchTerm.replace(/&/g, '%26');
            }
        }
        return params;// + '&appname=' + this.m_AppName;
    },

    //////////////////////////////////////////////////////////////////////////
    // gotoAndSearch - if term is HGVS expression, SNP, or position - find the
    // position requested. If it's search term - conduct the search.
    // parameters:
    // term - term to search or go to
    // options - Object with following possible fields
    //     view - if set, update this particular view, otherwise all views
    //     success - if set, report results through this callback,
    //               if not for position set it to the view, or all views if
    //               view is not specified. Set positional marker.
    //               For search, show search dialog.
    //         results - object of the following composition
    //             type - "hgvs",
    //             total - number of features found
    //             features - array of objects
    //                 from - first feature coordinate
    //                 to   - second feature coordinate
    //                 gi   - gi feture found on
    //                 label - feature label
    //         options - options Object passed to gotoAndSearch
    //     failure - if set, report errors through this callback
    //               if not, shows Message Boxes
    //         err_msg - error message text
    //         options - options Object passed to gotoAndSearch
    // returns array with deferred ids
    gotoAndSearch: function (term, options) {
        if (!term) {
            return;
        }
        term = Ext.util.Format.trim(term);
        if (!term.length) {
            return;
        }
        if (!options) {
            // Guard against empty options
            options = {};
        }
        // Add view and scope to callback so it works with parseAndGotoPosition
        if (!options.scope) options.scope = this;
        var hgvs = this.isHGVSExpression(term);
        var snip = term.match(/^([rs]s)([0-9]{3,})(?::.+)?$/);
        var vardb = term.match(/^([en]s)(td|v|sv)([0-9]+)(?::.+)?$/);
        if (hgvs) {
            return this.handleHGVS(term, options);
        } else if (snip) {
            return this.handleVariations(snip[1] + snip[2], options);
        } else if (vardb) {
            return this.handleVariations(vardb[1] + vardb[2] + vardb[3], options);
        } else if (this.isPosOrRange(term)) {
            return this.parseAndGotoPosition(term, options);
        } else {
            if (!options.success) {
                this.openSearchDlg(term, options.view);
                return [];
            } else {
                return this.startSearch(term, options);
            }
        }
    },


    isHGVSExpression: function (term) {
        return term.match(/^(?:([^:]+):)?([cgmnrp])\.(.+)$/);
    },

    //////////////////////////////////////////////////////////////////////////
    // setPositionalMarker - set a marker with a special property: initially locked
    // when unlocked it changes its title to default "Marker #"
    setPositionalMarker: function (pos_range, name) {
        if (!this.getMarkersInfo().findMarkerByName(name)
            && !(pos_range.length == 1 && this.getMarkersInfo().findMarkerByPos(pos_range[0]))
        ) {
            this.getMarkersInfo().addMarker(pos_range, name,
                SeqView.MarkerFlags.PositionTitle | SeqView.MarkerFlags.UserLock);
        }
    },




    //
    // term - position or range to parse and goto
    // options - Object with parameters:
    //     range_only - should we accept strictly range or position is OK, false
    //     view       - goto this position in a specific view, otherwise - globally
    //     success    - do not go to position, report successfully parsed position to this callback
    //         results in 'features' format
    //         options - options object for parseAndGotoPosition
    //     failure     - do not report error to user, use this callback
    //         error message
    //         options - options object for parseAndGotoPosition
    //     scope      - object to use as 'this' for callbacks
    parseAndGotoPosition: function (term, options) {
        term = term.replace(/[, ]/g, '');
        this.handlePos(term, {
            // pass parameters from our options
            range_only: options.range_only,
            view: options.view,
            nested_success: options.success,
            nested_failure: options.failure,
            // handlePos parameters
            ask_user: true,
            success: function (pos_range, options) {
                if (options.range_only && pos_range.length !== 2) {
                    options.failure.call(options.scope, "Range required", options);
                    return;
                }
                if (options.nested_callback) {
                    var pos = pos_range[0];
                    var pos_to = pos;
                    if (pos_range.length > 1) pos_to = pos_range[1];
                    var res = { type: "pos", total: 1,
                        features: [{ from: pos,
                            to: pos_to,
                            gi: options.app.m_ViewParams.gi,
                            label: term
                        }]
                    };
                    options.nested_callback.call(options.scope, res, options);
                } else {
                    if (pos_range.length == 1) {
                        options.app.setPositionalMarker(pos_range, term);
                    }
                    if (options.view) {
                        options.view.gotoPosRange(pos_range, true, { from_ui: true });
                    } else {
                        options.app.forEachView(function (v) {
                            if (v.isGraphic()) v.gotoPosRange(pos_range, true, { from_ui: true });
                        });
                    }
                }
            },
            failure: function (err_msg, options) {
                if (options.nested_failure) {
                    options.nested_failure.call(options.scope, err_msg, options);
                } else {
                    Ext.MessageBox.alert('Error', err_msg);
                }
            },
            scope: options.scope ? options.scope : this,
            app: this
        });
        return [];
    },

    checkSearchStatus: function (data, request_type, request_term, options) {
        var from_cgi = SeqView.decode(data);
        if (from_cgi.job_status) {
            var st = from_cgi.job_status;
            if (st == 'submitted' || st == 'running' || st == 'pending') {
                Ext.defer(SeqView.App.simpleAjaxRequest, 2000, this, [{
                    url: this.m_CGIs.FeatSearch + '?job_key=' + from_cgi.job_id,
                    context: this,
                    success: function (data, text, res) {
                        this.checkSearchStatus(data, request_type, request_term, options);
                    },
                    error: function (data, text, res) {
                        this.reportSearchFailure(request_type, request_term, options);
                    }
                }]);
                return;
            }
        }
        from_cgi.type = request_type;
        if (options.success)
            options.success.call(this, from_cgi);

    },
    reportSearchFailure: function (request_type, request_term, options) {
        var msg = "Server error in searching request " + request_term;
        if (options.failure)
            options.failure.call(this, msg);
    },

    startSearch: function (term, options) {
        var seq_type = 'nucleotide';
        if (this.m_ViewParams['acc_type'] == 'protein') {
            seq_type = 'protein';
        }
        var data_key = (this.m_App.m_Key && this.m_App.m_Key.length > 0) ? this.m_App.m_Key : '';
        var data_tracks = this.m_App.getActTracks();
        if (!data_tracks)
            data_tracks = '';

        var requests = [
            { url: this.m_CGIs.FeatSearch,
                params: { term: term, type: 'feature', id: this.GI, whole: 1, limit: 1000000, key: data_key, tracks: data_tracks }
            },
            { url: this.m_CGIs.FeatSearch,
                params: { term: term, type: 'component', id: this.GI, whole: 1, limit: 1000000, key: data_key }
            },
            { url: this.m_CGIs.FeatSearch,
                params: { term: term, type: seq_type, id: this.GI, whole: 1, limit: 1000000, key: data_key }
            }
        ];
        var res = [];

        for (var i = 0, l = requests.length; i < l; i++) {
            res.push((function (that, type, term) {
                return that.AjaxRequest({
                    url: requests[i].url,
                    data: requests[i].params,
                    context: that,
                    success: function (data, text, res) {
                        that.checkSearchStatus(data, type, term, options);
                    },
                    error: function (data, text, res) {
                        that.reportSearchFailure(type, term, options);
                    }
                });
            })(this, requests[i].params.type, requests[i].params.term));
        }
        return res;
    },

    /// MRU search pattern list support


    spsRecId: 100,
    searchPatternData: [],

    loadSearchPatternData: function () {
        var app = this;
        this.AjaxRequest({
            url: this.m_CGIs.SearchMru, context: this,
            data: {
                session: SeqView.SessionData.get("ncbi_sid") || '',
                seqid: this.GI,
                oper: 'load'
            },
            success: function (data, text, req) {
                if (data.mrulist) {
                    var mrulist = Ext.urlDecode(data.mrulist).list;
                    if (mrulist) {
                        Ext.each(Ext.util.JSON.decode(mrulist), function(pt) {app.addSearchPattern(pt, 'push');});
                        app.forEachView(function() {
                            var combo = this.m_View.down('#gotoBox');
                            if (combo) combo.setStore(app.searchPatternData);
                        });
                    }
                }
            },
            error: function () { }
        });

    },

    addSearchPattern: function(pattern, func) {
        if (typeof pattern != 'string') return;
        var maxNumOfSearchPatterns = 20;
        var spd = this.searchPatternData;
        var ix = spd.indexOf(pattern);
        if (ix >= 0) spd.splice(ix, 1); 
        spd[func || 'unshift'](pattern);
        spd.splice(maxNumOfSearchPatterns);
    },

    saveSearchPatternData: function () {
        var mrulist = Ext.urlEncode({list:  Ext.util.JSON.encode(this.searchPatternData)});
        this.AjaxRequest({
            url: this.m_CGIs.SearchMru, context: this,
            data: {
                session: SeqView.SessionData.get("ncbi_sid") || '',
                seqid: this.GI,
                oper: 'save',
                mrulist: mrulist
            },
            success: function (data, text, req) { },
            error: function (data, text, req) { }
        });
    },

    showSearchParamsDlg: function (view) {
        var combo = Ext.create('Ext.form.field.ComboBox', {
            listeners: {
                specialkey: function(f,e){
                    if (e.getKey() == e.ENTER){
                       this.searchRequest();
                    }
                }
            },
            typeAhead: true,
            store: this.searchPatternData,
            
           searchRequest: function() { 
               var term = combo.getValue();
               if (!term) return;
               dlg.close();
               view.gotoAndSearch(term);
           }
        });
        var dlg = new Ext.Window({
            title: 'Find on Sequence',
            app: this,
            modal: true,
            layout: 'form',
            width: 420,
            padding: 3,
            buttonAlign: 'center',
            buttons: [
                {text: 'OK', handler: combo.searchRequest},
                {text: 'Cancel', handler: function(){ dlg.close(); }}]
        });
        dlg.add({xtype: 'displayfield', bodyStyle: {background: 'inherit'}, border: false,
            value: 'Enter sequence position or range (possible range formats are '
            + '10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m)<br/><br/>'
            + 'Or enter name of feature, component, HGVS, SNP rs id, '
            + 'or sequence (nucleotide regexp with IUPAC equivalents, PROSITE patterns):<br/>'})

        dlg.add(combo);
        dlg.show();
    },

    showSearchDlg: function () {
        var str = this.m_View.down('#gotoBox').getValue();
        openSearchDlg(str);
    },

    openSearchDlg: function (seek, view) {
        var str = seek;
        if (str.length == 0) {
            str = '*';
        }

        if (this.m_SearchDlg) {
            this.m_SearchDlg.close();
            this.m_SearchDlg = null;
        }
        var fixDelay = function(ajax, res, opt) {
            var from_cgi = Ext.JSON.decode(res.responseText);
            opt.callback = opt.callback || opt.callbackBAK;
            if (from_cgi.job_status && from_cgi.job_id) {
                opt.callbackBAK = opt.callback || opt.callbackBAK;
                delete opt.callback;
                opt.params = {job_key: from_cgi.job_id};
                Ext.defer(Ext.Ajax.request, 2000, this, [opt]);
           } 
           return true;
        }
                
        if (!Ext.Ajax.hasListeners.requestcomplete) {
            Ext.Ajax.on('requestcomplete', fixDelay) 
        }
        
        var searchTab = function (title, view, str, tabpanel) {
            this.m_View = view;
            var params = { term: str, id: this.m_View.m_App.GI, from: this.m_View.m_FromSeq,
                to: (this.m_View.m_FromSeq + this.m_View.m_LenSeq - 1),
                limit: 50, whole: 1, appname: this.m_View.m_App.m_AppName
            };
            if (this.m_View.m_App.m_Key && this.m_View.m_App.m_Key.length > 0) {
                params.key = this.m_View.m_App.m_Key;
            }
            switch (title) {
                case "Components":
                    params.type = 'component';
                    break;
                case "Sequence":
                    if (view.m_App.m_ViewParams['acc_type'] == 'protein') {
                        params.type = 'protein';
                    } else {
                        params.type = 'nucleotide';
                    }
                    break;
                case "Features":
                    params.type = 'feature';
                    var tracks = this.m_View.m_App.getActTracks();
                    if (tracks)
                        params.tracks = tracks;
                    break;
            }
            var proxy = {
                type: 'ajax',
                url: view.m_App.m_CGIs.FeatSearch, extraParams: params,
                reader: { rootProperty: 'features', totalProperty: 'total' },
                simpleSortMode: true
            };
            var fields = [
                'label',
                {name: 'from', type: 'int'},
                {name: 'to', type: 'int'},
                {name: 'strand', mapping: function(data)
                    {return data.strand == '-' ? 'Negative' : 'Positive';}}
            ];
            


            var store = Ext.create('Ext.data.Store', {
                model: Ext.define(null, {
                    extend: 'Ext.data.Model',
                    fields: fields,
                    proxy: proxy }),
                pageSize: 50,
                remoteSort: true,
                autoLoad: true,
                job_id: null,
                sorters: [{ property: 'label', direction: 'ASC' }],
                listeners: {
                    load: function (obj, data) {
                        if (this.panel.isDestroyed) return;
                        this.panel.setIconCls(data && data.length > 0 ? 'xsv-search-results' : '');
                }}
            });
            var panel = new Ext.grid.GridPanel({
                title: title,
                store: store,
                iconCls: 'xsv-search-loading',
                border: false,
//                loadMask: true,
                stripeRows: true,
                columns:[
                    { text: "Label", flex: true, wi_dth: 110, sortable: true, dataIndex: 'label' },
                    { text: "From", width: 90, sortable: true, dataIndex: 'from', renderer: function (val) { return val + 1; } },
                    { text: "To", width: 90, sortable: true, dataIndex: 'to', renderer: function (val) { return val + 1; } },
                    { text: "Strand", width: 90, sortable: true, dataIndex: 'strand'} 
                ],
                viewConfig: {forceFit:true, deferEmptyText:false, emptyText:'<div align="center">No Search Results To Display</div>'},
                listeners: {
                    'rowdblclick': function (table, rec, tr, rowIndex, e) {
                        var view = this.m_View;
                        if (table.ownerGrid.title == 'Component') {
                            var tracks_changed = view.m_App.showTracks({ key: 'component_track' }, true, true);
                            tracks_changed |= view.m_App.showTracks({ key: 'scaffold_track' }, true, true);
                            if (tracks_changed) view.m_App.fireEvent('configuration_changed', view.m_App);
                        }
                        var from = rec.get('from');
                        var to = rec.get('to');
                        var view_len = to - from + 1;

                        view.m_SelectedSig = rec.get('object_id') || null;
                        if (view.m_SelectedSig)
                            view.removeRangeSelection(true);
                        else
                            view.m_SelectedRangeSet =[[from, to]];
                        view.pingClick('2-3-' + table.ownerGrid.title.charAt(0));
                        view.startImageLoading(from - view_len * .15, view_len * 1.3);
                        view.m_App.fireEvent('selection_changed', view); 
                    },
                    scope: this
                },
                bbar: new Ext.PagingToolbar({
                    store: store,
                    displayInfo: true, displayMsg: 'Displaying Search Results {0} - {1} of {2}', emptyMsg: "No Search Results To Display"
                })

            });
            store.panel = panel;
            panel.on('beforedestroy', function () {
                if (this.store.job_id != null) {
                    SeqView.App.simpleAjaxRequest({ url: this.store.url + '?Cancel=true&job_key=' + this.store.job_id });
                }
            });
            return panel;
        };

        // This is a hack. we need to find a first graphic view in order to get the search range
        // ideally, the search function should be called with the view parameter
        this.forEachView(function (v) { if (!view && v.isGraphic()) view = v; });
        if (view) {
            var tabPanel = new Ext.TabPanel({
                activeTab: 0,
                listeners: {
                    'tabchange': function (tabpanel, tab) {
                        tabpanel.down('#reconf_btn').setVisible(tab.title == 'Tracks');
                    }
                },
                buttons: [
                    {text: 'Re-configure', hidden: true, scope: this, itemId: 'reconf_btn',
                        handler: function () {
                            this.pingClick('2-3-T');
                            if (this.m_SearchDlg.getActiveTab().syncModifiedTracks()) {
                                var categories = SeqView.TM.processTracksInfo(this.m_Config.TrackConfig);
                                SeqView.TM.Common.updateSeqViewApp(categories, this);
                            }
                            this.m_SearchDlg.close(); this.m_SearchDlg = null;
                        }
                    },
                    {text: 'Close', scope: this, handler: function () { this.m_SearchDlg.close(); this.m_SearchDlg = null; } }
                ]
            });
            tabPanel.add(new searchTab("Features", view, str, tabPanel));
            tabPanel.add(new searchTab("Components", view, str, tabPanel));
            tabPanel.add(new searchTab("Sequence", view, str, tabPanel));
            tabPanel.add(SeqView.TM.createSearchTracksTab("Tracks", this.m_Config.TrackConfig, str));
            this.m_SearchDlg = new Ext.Window({
                app: this,
                layout: 'fit',
                title: 'Search Results',
                minWidth: 480, width: 650, height: 350,
                constrain: true, collapsible: true,
                cls: 'SeqViewerApp',
                closeAction: 'destroy',
                plain: true,
                items: [tabPanel],
                getActiveTab: function ()
                { return tabPanel.activeTab; }
            });
            var app = this;
            this.m_SearchDlg.on('close', function () {
                app.resizeIFrame(); app.m_DialogShown = false;
            });

            app.resizeIFrame(400);
            app.m_DialogShown = true;
            this.m_SearchDlg.show();
        }
    },

    showFeedbackDlg: function (logarea, fbText, fbURL) {
        this.pingClick(logarea);
        var app = this;
        if (typeof fbText != 'string') fbText = null;
        var fbCallback = function (linkUrl) {
            app.resizeIFrame(400);
            var feedbackType = ['Suggestion', 'Bug Report', 'Other', 'Initial upload error'];
            var fbType = (fbText) ? 3 : 0;
            if (!fbType) feedbackType.pop();
            var feedbackDlg = new Ext.Window({
                app: this,
                layout: 'fit', modal: true,
                title: 'NCBI Graphical Sequence Viewer feedback',
                width: 600, height: 360,
                minWidth: 400, minHeight: 260,
                resizable: true,
                cls: 'SeqViewerApp',
                closeAction: 'destroy',
                items: [{
                    xtype: 'form',
                    itemId: 'fbForm',
                    bodyStyle: {background: 'inherit'},
                    labelWidth: 140,
                    frame: false,
                    labelAlign: 'right',
                    items: [
                         {xtype: 'combo', triggerAction: 'all', fieldLabel: 'Feedback Type', mode: 'local', name: 'feedback-type',
                              store: feedbackType, allowBlank: false, editable: false, value: feedbackType[fbType], disabled: (fbType > 0)},
                         {xtype: 'textfield', fieldLabel: 'EMail (Optional)', vtype: 'email', allowBlank: true, anchor: '100%', name: 'feedback-email' },
                         {xtype: 'textarea', fieldLabel: '*Feedback', allowBlank: false, anchor: '100% -50', name: 'feedback-text', value: fbText }
                      ]
                }],
                buttons: [
                  { text: 'Send', handler: function () {
                      var form = feedbackDlg.down('#fbForm').getForm();
                      if (form.isValid()) {
                          var cfg_data = form.getValues();
                          cfg_data['feedback-browser'] = Ext.browser.name + ' ver. ' + Ext.browser.version.version;
                          cfg_data['feedback-os'] = Ext.os.name;
                          NCBIGBUtils.makeTinyURL(linkUrl, function (res) {
                              cfg_data['feedback-url'] = linkUrl + '\n' + (res.id || '');
                              app.AjaxRequest({ url: app.m_CGIs.Feedback, data: cfg_data, context: feedbackDlg,
                                  success: function (data) {
                                      this.close();
                                      Ext.MessageBox.show({ title: 'Feedback', msg: 'Thank you! We appreciate your feedback.',
                                          buttons: Ext.MessageBox.OK, icon: Ext.MessageBox.INFO
                                      });
                                  },
                                  error: function (data, txt, res) {
                                      console.log('Failed to send feedback: ' + txt);
                                  }
                              });
                          });
                      }
                  }
                  },
                  { text: 'Cancel', handler: function () { feedbackDlg.close(); } }
               ]
            });
            feedbackDlg.on('close', function () { app.resizeIFrame(); app.m_DialogShown = false; });

            app.m_DialogShown = true;
            if (app.m_iFrame) Ext.defer(feedbackDlg.show, 500, feedbackDlg); else feedbackDlg.show();
        }
        if (fbText) fbCallback(fbURL);
        else this.getLinkToThisPageURL(null, null, fbCallback);
    },

    findSelectionToolTip: function(key) {
        return this.m_SelectionTTMap.get(key);
    },

    saveSelectionToolTip: function(key, value) {
        this.m_SelectionTTMap.set(key, value);
    }

});

//////////////////////////////////////////////////////////////////////////
// SeqView.PrinterFriendlyPage


SeqView.PrinterFriendlyPage = function(app, printdlg) {
    this.values = {orientation: 1, sequence: true, title: true}
    var img_width = this.values['orientation']==2 ? 900 : 650;
    var marker_num = 0;
    this.generator=window.open('','name','width='+(img_width+50)+',height=700,toolbar=0,location=0,status=0,menubar=1,resizable=1,scrollbars=1');
    this.generator.document.write('<html><head><title>'+app.GI+'</title>');
    this.generator.document.write('<link rel="stylesheet" href="'+app.m_CGIs.prefix+'css/style.css">');
    var callFixSeqTrans = false;
    if (this.values['sequence'] && app.m_TextView) {
        if (app.m_TextView.m_HideTrans) {
            callFixSeqTrans = true;
            this.generator.document.write('<script>function fixSeqTrans() {');
            this.generator.document.write('var arr = document.getElementsByTagName("div");');
            this.generator.document.write('for (var i = 0; i < arr.length; i++) {');
            if (app.m_TextView.m_HideTrans == 1) {
                this.generator.document.write('if (arr[i].className.indexOf("seqtrans_trans") == 0) {');
            } else if (app.m_TextView.m_HideTrans == 2) {
                this.generator.document.write('if (arr[i].className.indexOf("seqtrans_prot") == 0) {');
            } else if (app.m_TextView.m_HideTrans == 3) {
                this.generator.document.write('if (arr[i].className.indexOf("seqtrans") == 0) {');
            }
            this.generator.document.write('arr[i].style.display = "none";');
            this.generator.document.write('}}}</script>');
        }
    }
    this.generator.document.write('</head><body onLoad="' + (callFixSeqTrans ? "fixSeqTrans();" : "") + 'window.print()">');
    if (this.values['title'])  {
        this.generator.document.write('<p><b>'+app.m_ViewParams['title']+'</b></p>'); // add title
        this.generator.document.write('<p>'+app.m_ViewParams['id_full']+'</p>'); // add full id
    }

    this.ondataloadcalled = false;
    this.running_requests = 0;
    this.overview_data = null;

    var view_data_idx = 0;
    this.views_data = [];

    this.text_data = null;
    if (this.values['sequence'] && app.m_TextView) { // checkbox checked
        this.col_tmp = app.m_TextView.m_SequenceCols;
        app.m_TextView.m_SequenceCols = this.values['orientation'] == 2 ? (Ext.isWindows ? 100: 110) : (Ext.isWindows ? 70 :80);
        if( app.m_TextView.m_Flip ){
            var numLines = app.m_TextView.m_LenSeq / app.m_TextView.m_SequenceCols;
            numLines = Math.ceil( numLines );
            app.m_TextView.m_LenSeq = numLines * app.m_TextView.m_SequenceCols;
        }
        var the_url = app.m_TextView.getSeqTextURL();
        
        this.running_requests++;
        app.AjaxRequest({url: the_url, context: this,
            success:function(data, text, res) {
                this.text_data = SeqView.decode(data);
                this.running_requests--;
                this.onDataLoaded();
            }
       });
    }

    if ( this.ondataloadcalled === false && this.running_requests == 0) {
        // only sequence title is printed
        this.generator.document.write('</body></html>');
        this.generator.document.close();
    }

    this.onDataLoaded = function() {
        this.ondataloadcalled = true;

        if (this.running_requests != 0)
            return;
        if (this.overview_data) {
            // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
            // This is a way to provide reliable URL resolution for embedding. SV-1760
            var img_url = this.overview_data.img_url;
            if (img_url && img_url.charAt(0) == '?') {
                img_url = app.m_CGIs.NetCache + img_url;
            }
            this.generator.document.write('<img style="border:1px solid;" src="' + img_url + '">');
        }
        var print_tpl_str = '<div style="position:absolute;top:2px;width:16px;height:16px;left:{left}px;"><div class="marker_label" style="{font_size}border:1px {color} solid;color:black;">{trimmed_label}</div><div class="marker_line" style="border-right:1px solid {color}; height:{height}px;"></div></div>'
        Ext.each(this.views_data, function(view_data) {
            var view = view_data.view;
            var from_cgi = view_data.data;
            this.generator.document.write('<br><br>'+view.m_View.title + '<br>');
            var html = '<div style="position:relative;">';
            var margin = 0;
            if (app.m_MarkersInfo) {
                app.m_MarkersInfo.forEachMarker( function(the_m) {
                    if ( the_m.seq_pos>view.m_VisFromSeq && the_m.seq_pos<(view.m_VisFromSeq+view.m_VisLenSeq-1)) {
                        var bpp = from_cgi.len / img_width;
                        var marker_left = 0;
                        if (view.getFlip())
                            marker_left = (from_cgi.len - (the_m.seq_pos - from_cgi.from) -1 + 0.5) / bpp - 8;
                        else
                            marker_left = (the_m.seq_pos + 0.5 - from_cgi.from) / bpp - 8;

                        var options = {color:the_m.color, left:marker_left, height:from_cgi.h + 4,
                                    num:marker_num, trimmed_label:the_m.marker_name.trimToPix(100),
                                    font_size:'font-size:0.95em;'
                        };
                        marker_num += 1;

                        var print_tpl = new Ext.Template(print_tpl_str);
                        var marker_html = print_tpl.applyTemplate(options);
                        html += marker_html;
                    }
                });
                margin = 21;
            }
            // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
            // This is a way to provide reliable URL resolution for embedding. SV-1760
            var img_url = from_cgi.img_url;
            if (img_url && img_url.charAt(0) == '?') {
                img_url = app.m_CGIs.NetCache + img_url;
            }
            html += '<img style="border: 1px solid;margin-top:'+margin+'px;" src="' + img_url + '"></div>';
            this.generator.document.write(html);


        },this);

        if (this.text_data) {
            this.generator.document.write('<br><br>Sequence View<br>');
            var top_nums = app.m_TextView.genTopNumbers();
            var html = SeqView.TextView.getViewTmpl().apply({border:1, width:img_width, left_nums:this.text_data.starts, top_nums:top_nums, id: app.m_Idx, sequence:this.text_data.sequence});
            this.generator.document.write(html);
            app.m_TextView.m_SequenceCols = this.col_tmp; // restore back
        }

        this.generator.document.write('</body></html>');
        this.generator.document.close();

    }

};
function globStringToRegexStr(str) {
    return preg_quote(str).replace(/\\\*/g, '.*').replace(/\\\?/g, '.');
}

function globStringToRegex(str) {
    return new RegExp(preg_quote(str).replace(/\\\*/g, '.*').replace(/\\\?/g, '.'), 'i');
}
function preg_quote (str, delimiter) {
    // http://kevin.vanzonneveld.net
    // +   original by: booeyOH
    // +   improved by: Ates Goral (http://magnetiq.com)
    // +   improved by: Kevin van Zonneveld (http://kevin.vanzonneveld.net)
    // +   bugfixed by: Onno Marsman
    // +   improved by: Brett Zamir (http://brett-zamir.me)
    // *     example 1: preg_quote("$40");
    // *     returns 1: '\$40'
    // *     example 2: preg_quote("*RRRING* Hello?");
    // *     returns 2: '\*RRRING\* Hello\?'
    // *     example 3: preg_quote("\\.+*?[^]$(){}=!<>|:");
    // *     returns 3: '\\\.\+\*\?\[\^\]\$\(\)\{\}\=\!\<\>\|\:'
    return (str + '').replace(new RegExp('[.\\\\+*?\\[\\^\\]$(){}=!<>|:\\' + (delimiter || '') + '-]', 'g'), '\\$&');
}
}
