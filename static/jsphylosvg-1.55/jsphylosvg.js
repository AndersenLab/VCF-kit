Smits = {};Smits.Common = {
	nodeIdIncrement : 0,
	activeNode: 0,
	
	/* Rounds float to a defined number of decimal places */
	roundFloat : function(num, digits){
		var i = 0, 
			dec = 1;
		while(i < digits){
			dec *= 10;
			i++;
		}
		return Math.round(num*dec)/dec; 
	},
	
	/* Copies properties from one object to another */
	apply : function(obj, extObj){
		if (obj && typeof extObj == 'object') {
			for (var key in extObj) {
				obj[key] = extObj[key];
			}
		}
		return obj;	
	},
	
	addRaphEventHandler : function(el, eventType, fn, paramsObj){
		try{
			el[eventType](function(fn, paramsObj){
				return function(e,o){
					var params = paramsObj;
					params.e = e;
					fn(params);
				};
			}(fn, paramsObj));
		} catch (err){}	
	},
	
	isInteger : function(s) {
		return !isNaN(parseInt(s));
	},

	isXMLSerializerAvailable : function(){
		if (typeof(XMLSerializer) == "function"){
			return true;
		} else {
			return false;
		}
	},
	
	createSvgEl : function (el, attr) {
		el = document.createElementNS("http://www.w3.org/2000/svg", el);            
		if (attr) {
			for (var key in attr) {
				if (attr.hasOwnProperty(key)) {
					el.setAttribute(key, String(attr[key]));
				}
			}
		}	
		return el;	
	},
	
	createGradientEl : function(name, obj, coords){
		if(obj.type != "radialGradient") return false;
		
		var radialEl = Smits.Common.createSvgEl("radialGradient", {
			id: name, 
			gradientUnits:"userSpaceOnUse", 
			cx: coords[0], 
			cy: coords[1], 
			r: coords[2], 
			fx: coords[0], 
			fy: coords[1]
		});

		if(obj.stop){
			var stop = obj.stop;
			for(var i = 0; i < stop.length; i++){
				var stopObj = stop[i];
				if(stopObj['@attributes']){
					radialEl.appendChild(Smits.Common.createSvgEl("stop", stopObj['@attributes']));
				} else {
					if(stopObj['_attributes']) delete stopObj['_attributes'];
					if(stopObj['_children']) delete stopObj['_children'];
					if(stopObj['__proto__']) delete stopObj['__proto__'];
					radialEl.appendChild(Smits.Common.createSvgEl("stop", stopObj));
				}
			}
		}
		
		return radialEl;
	},
	
	setCssStyle : function(selector, rule) {
		var stylesheet = document.styleSheets[0];
		if( stylesheet.addRule ){
			stylesheet.addRule(selector, rule);
		} else if( stylesheet.insertRule ){
			stylesheet.insertRule(selector + ' { ' + rule + ' }', stylesheet.cssRules.length);
		}
	}

};Smits.PhyloCanvas = function(){
	var phylogram,
		divId,
		newickObject,
		svg,
		dataObject;

	return function(inputFormat, sDivId, canvasWidth, canvasHeight, type){
		/* Privileged Methods */
		this.getNewickObject = function(){
			return newickObject;
		};
		this.clear = function(){
		
		};
		this.scale = function(multiplier){
			svg.svg.scale(multiplier);
		};
		this.getSvg = function(){
			return svg;
		};
		this.getPhylogram = function(){
			return phylogram;
		};
		this.getSvgSource = function(){
			if(Raphael.svg && Smits.Common.isXMLSerializerAvailable()){
				var serialize = new XMLSerializer();
				return serialize.serializeToString(svg.svg.canvas);
			} else {
				return false;
			}
		}
	
		/* CONSTRUCTOR */

		// Process dataset -- assume newick format, else needs to provide format
		if(typeof inputFormat === "object"){
			if(inputFormat.xml){	// default xml format is phyloXML
				if(!inputFormat.fileSource){
					var xj = XMLObjectifier.textToXML(inputFormat.xml); 			// assume we need to clean it up
				} else {
					var xj = inputFormat.xml;
				}
				xj = XMLObjectifier.xmlToJSON(xj);
				dataObject = new Smits.PhyloCanvas.PhyloxmlParse(xj);
			} else if(inputFormat.phyloxml){
				if(!inputFormat.fileSource){
					var xj = XMLObjectifier.textToXML(inputFormat.phyloxml); 			// assume we need to clean it up
				} else {
					var xj = inputFormat.phyloxml;
				}
				xj = XMLObjectifier.xmlToJSON(xj);
				dataObject = new Smits.PhyloCanvas.PhyloxmlParse(xj);
			} else if(inputFormat.nexml){
				if(!inputFormat.fileSource){
					var xj = XMLObjectifier.textToXML(inputFormat.nexml); 			// assume we need to clean it up
				} else {
					var xj = inputFormat.nexml;
				}
				xj = XMLObjectifier.xmlToJSON(xj);
				dataObject = new Smits.PhyloCanvas.NexmlParse(xj, inputFormat)
			} else if(inputFormat.json){
				dataObject = new Smits.PhyloCanvas.PhyloxmlParse(inputFormat.json);
			} else if(inputFormat.newick){
				dataObject = new Smits.PhyloCanvas.NewickParse(inputFormat.newick);
			} else if(inputFormat.nexmlJson){
				dataObject = new Smits.PhyloCanvas.NexmlJsonParse(inputFormat);				
			} else {
				alert('Please set the format of input data');
			}
		} else {
			dataObject = new Smits.PhyloCanvas.NewickParse(inputFormat);
		}

		divId = sDivId;
		svg = new Smits.PhyloCanvas.Render.SVG( divId, canvasWidth, canvasHeight );
		
			/* FACTORY */
		if(type == "circular"){
			phylogram = new Smits.PhyloCanvas.Render.CircularPhylogram(
				svg, 
				dataObject
			);		
		} else {
			phylogram = new Smits.PhyloCanvas.Render.Phylogram(
				svg,
				dataObject
			);			
		}		
		
	}
	
}();

Smits.PhyloCanvas.prototype = {
};Smits.PhyloCanvas.Node = function(){
	/**
	* Node Class
	* Allows objects to be traversed across children
	*
	*/
	return function(o, parentInstance){
		// initiate object
		this.id = Smits.Common.nodeIdIncrement += 1;
		this.level = 0;
		this.len = 0;
		this.newickLen = 0;
		this.name = '';
		this.type = '';
		this.chart = {};
		this.img = [];
		
		if(o) Smits.Common.apply(this, o);

		/* Cache Calculations */
		this._countAllChildren = false;
		this._countImmediateChildren = false;
		this._midBranchPosition = false;
		
		this.children = new Array();
		
		if(parentInstance){
			parentInstance.children.push(this); 
		}
	}
}();


Smits.PhyloCanvas.Node.prototype = {
	
	getCountAllChildren : function(){
		if( this._countAllChildren !== false ) return this._countAllChildren;
		var nodeCount = 0;

		for (var key in this.children) {
			if(Smits.Common.isInteger(key)){
				var child = this.children[key];
				if(child.children && child.children.length > 0){
					nodeCount += child.getCountAllChildren();
				} else {
					nodeCount ++;
				}			
			}
		}
		this._countAllChildren = nodeCount;
		return nodeCount;
	},
	
	getCountImmediateChildren : function(){
		if( this._countImmediateChildren !== false ) return this._countImmediateChildren;
		var nodeCount = 0;

		for (var key in this.children) {
			var child = this.children[key];
			nodeCount += child.length;
		}
		this._countImmediateChildren = nodeCount;
		return nodeCount;
	},
	
	getMidbranchPosition : function(firstBranch){
		if( this._midBranchPosition !== false ) return this._midBranchPosition;
		var y = [0,0];  // bounds
		
		for (var i = 0; i < this.children.length; i++) {
			var child = this.children[i];
			if(child.children && child.children.length > 0){
				if(i == 0 && firstBranch){
					y[0] = child.getMidbranchPosition(true);				
					y[1] += child.getCountAllChildren() - 1;	
				} else if(i == 0){
					y[0] = child.getMidbranchPosition();				
					y[1] += child.getCountAllChildren();	
				} else if (i == this.children.length - 1){
					y[1] += child.getMidbranchPosition();
				} else {
					y[1] += child.getCountAllChildren();				
				}
			} else {
				if(i == 0 && firstBranch){
					y[0] = 0;
				} else if(i == 0){
					y[0] = 1;
					y[1] += 1;	
				} else if (i == this.children.length - 1){
					y[1] += 1;
				} else {
					y[1] += 1;
				}
			}
		}
		
		this._midBranchPosition = y[1] >= y[0] ? ((y[1] + y[0]) / 2) : y[0];
		return this._midBranchPosition;
	}
	
};Smits.PhyloCanvas.NewickParse = function(){

	var text,
	ch,
	pos,
	mLevel = 0,
	mNewickLen = 0,
	root,
	validate,
		
	object = function (parentNode) {
		var node  = new Smits.PhyloCanvas.Node();
		
		while (ch !== ')' && ch !== ',') {
			if (ch === ':'){
				next();
				node.len = Smits.Common.roundFloat(string(), 4);			// round to 4 decimal places
				if(node.len == 0){
					node.len = 0.0001;
				}
			} else if (ch === "'" || ch === '"'){ 
				node.type = "label";
				node.name = quotedString(ch);
			} else {
				node.type = "label";
				node.name = string();
			}
		}
		node.level = parentNode.level + 1;
		return node;
	},
	
	objectIterate = function(parentNode){
		var node = new Smits.PhyloCanvas.Node();
		if(parentNode){
			node.level = parentNode.level + 1;
		}
		
		while( ch !== ')' ){
			next();
			if( ch === '(' ) {
				node.children.push(objectIterate(node));
			} else {
				node.children.push(object(node));
			}
		}
		
		next();
		if(ch !== ':' && ch !== ')' && ch !== ',' && ch !== ';'){
			node.type = "label";
			node.name = string();
		}
		if(ch === ':'){
			next();
			node.len = Smits.Common.roundFloat(string(), 4);
			if(node.len == 0){
				node.len = 0.0001;
			}
			node.type = "stem";

		}
		return node;		
	},
	
	string = function(){
		var string = '';
		
		while (ch !== ':' && ch !== ')' && ch !== ',' && ch !== ';'){
			string += ch;
			next();
		}
		return string;
	},

	quotedString = function(quoteType){
		var string = '';
		
		while (ch !== quoteType){
			string += ch;
			next();
		}
		return string;
	},	
	
	next = function() {
		ch = text.charAt(pos);
		pos += 1;
		return ch;
	},
	
	recursiveProcessRoot = function(node, parentNode){
		
		if(node.children && node.children.length){
			for( var i = 0; i < node.children.length; i++ ){
				var child = node.children[i];
				if(child.len === 0) {	// Dendogram
					child.len = 1;	
				}
				child.newickLen = Smits.Common.roundFloat(child.len + node.newickLen, 4);
				if(child.level > mLevel) mLevel = child.level;
				if(child.newickLen > mNewickLen) mNewickLen = child.newickLen;
				if(child.children.length > 0){
					recursiveProcessRoot(child, node); 
				}				
			}
		}
		return node;
	};

	return function(parseText){
		/* Privileged Methods */
		this.getRoot = function(){
			return root;
		};
		this.getLevels = function(){
			return mLevel;
		};
		this.getNewickLen = function(){
			return mNewickLen;
		};		
		this.getValidate = function(){
			return validate;
		};		
		
		
		/* CONSTRUCTOR */	
		mLevel = 0;
		mNewickLen = 0;
		
		text = parseText;
		pos = 0;
		
		next();
		root = objectIterate();
		root = recursiveProcessRoot(root);
	}

}();

Smits.PhyloCanvas.NewickParse.prototype = {

};Smits.PhyloCanvas.PhyloxmlParse = function(){

	var mLevel = 0,
	mNewickLen = 0,
	root,
	validate,
		
	recursiveParse = function(clade, parentNode){
		var node = new Smits.PhyloCanvas.Node();
		if(parentNode){
			node.level = parentNode.level + 1;
		}
		
		if(clade.clade && clade.clade.length){
			for(var i = 0; i < clade.clade.length; i++){
				var thisClade = clade.clade[i];
				node.children.push(recursiveParse(thisClade, node));
			}
		}
		if(clade.branch_length){	// Branches can be attributes or own element
			if(typeof clade.branch_length === 'object'){
				clade.branch_length = clade.branch_length[0].Text;
			}

			node.len = Smits.Common.roundFloat(clade.branch_length, 4);			// round to 4 decimal places
			if(node.len == 0){
				node.len = 0.0001;
			}			
		}
		if(clade.name){
			node.type = 'label';
			node.name = clade.name[0].Text;
			if(clade.name[0] && clade.name[0].style){
				node.style = clade.name[0].style;
			}
			if(clade.name[0] && clade.name[0].bgStyle){
				node.bgStyle = clade.name[0].bgStyle;
			}			
		} else if(clade.confidence){
			node.name = clade.confidence[0].Text;
		}

		/* Collect further info that might be used as a label */
		if (clade.sequence && clade.sequence[0] && clade.sequence[0].name && clade.sequence[0].name[0] && clade.sequence[0].name[0].Text){
			node.sequenceName = clade.sequence[0].name[0].Text;
		}
		if (clade.taxonomy && clade.taxonomy[0]){
			if(clade.taxonomy[0].scientific_name && clade.taxonomy[0].scientific_name[0] && clade.taxonomy[0].scientific_name[0].Text){
				node.taxonomyScientificName = clade.taxonomy[0].scientific_name[0].Text;
			}
			if (clade.taxonomy[0].common_name  && clade.taxonomy[0].common_name[0] && clade.taxonomy[0].common_name[0].Text){
				node.taxonomyCommonName = clade.taxonomy[0].common_name[0].Text;
			}
		}
		if (clade.sequence && clade.sequence[0] && clade.sequence[0].accession && clade.sequence[0].accession[0] && clade.sequence[0].accession[0].Text){
			node.sequenceAccession = clade.sequence[0].accession[0].Text;
		}
		if (clade.point ){
			node.LatLong = [clade.point[0].lat[0].Text, clade.point[0]['long'][0].Text];
		}		

		
		/* Prioritization of Label */
		if(!node.name){
			if(node.sequenceName){
				node.name = node.sequenceName;
			} else if (node.taxonomyScientificName){
				node.name = node.taxonomyScientificName;
			} else if (node.taxonomyCommonName){
				node.name = node.taxonomyCommonName;
			} else if (node.sequenceAccession){
				node.name = node.sequenceAccession;
			}
			if(node.name){	// if name is now set, type is 'label'
				node.type = 'label'; 
			}
		}
		
		if(clade.annotation){
			if(clade.annotation[0] && clade.annotation[0].desc && clade.annotation[0].desc[0] && clade.annotation[0].desc[0].Text){
				node.description = clade.annotation[0].desc[0].Text;
			}
			if(clade.annotation[0] && clade.annotation[0].uri && clade.annotation[0].uri[0] && clade.annotation[0].uri[0].Text){
				node.uri = clade.annotation[0].uri[0].Text;
			}			
			if(clade.annotation[0] && clade.annotation[0].img){
				for(var i in clade.annotation[0].img){
					if(Smits.Common.isInteger(i)){
						node.img[i] = clade.annotation[0].img[i].Text;
					}
				}
			}
		}
		if(clade.chart){
			if(clade.chart[0]){
				for(var i in clade.chart[0]){
					if(i != 'Text' && i != '_children'){
					node.chart[i] = clade.chart[0][i][0].Text;
					}
				}
			}
			
		}
		
		// Validation
		if(node && node.level > 1){
			if(!node.len){
				validate = 'Error. Please include Branch Lengths - we only draw rooted phylogenetic trees.';
			}
		}
			
		return node;
	},
	
	recursiveProcessRoot = function(node, parentNode){
		if(node.children && node.children.length){
			for( var i = 0; i < node.children.length; i++){
				var child = node.children[i];
				child.newickLen = Math.round( (child.len + node.newickLen) *10000)/10000;
				if(child.level > mLevel) mLevel = child.level;
				if(child.newickLen > mNewickLen) mNewickLen = child.newickLen;
				if(child.children.length > 0){
					recursiveProcessRoot(child, node); 
				}				
			}
		}
		return node;
	},
	
	recursiveProcessParameters = function(parametersEl, treeType){
		for (var i in parametersEl){
			if(i != '_children' && i != 'Text'){
				if(i == 'rectangular' || i == 'circular'){
					recursiveProcessParameters(parametersEl[i][0], i);
				} else {
					if(!Smits.PhyloCanvas.Render.Parameters[i]) {  Smits.PhyloCanvas.Render.Parameters[i] = {}; };
					Smits.PhyloCanvas.Render.Parameters.set(i, parametersEl[i][0].Text, treeType);
				}
			}
		}
		return;
	};

	return function(jsonString){
		/* Privileged Methods */
		this.getRoot = function(){
			return root;
		};
		this.getLevels = function(){
			return mLevel;
		};
		this.getNewickLen = function(){
			return mNewickLen;
		};		
		this.getValidate = function(){
			return validate;
		};
		
		
		/* CONSTRUCTOR */	
		if(jsonString.phylogeny && jsonString.phylogeny[0] && jsonString.phylogeny[0].clade){
			root = recursiveParse(jsonString.phylogeny[0].clade[0]);
		}
		
		if(jsonString.phylogeny && jsonString.phylogeny[0] && jsonString.phylogeny[0].render && jsonString.phylogeny[0].render[0]){
			var render = jsonString.phylogeny[0].render[0];
			
			// Custom Styles
			if(render && render.styles){
				var styles = render.styles[0];
				for (var i in styles){
					if(i != '_children' && i != 'Text'){
						if(styles[i][0]['type'] && styles[i][0]['type'] == "radialGradient" && Raphael.svg){
							// radialGradient only supported by SVG
							styles[i][0]['name'] = i;
							Smits.PhyloCanvas.Render.Style[i] = styles[i][0];
							if(!Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList']) { Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'] = [] };
							Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'].push(i); 
						} else {
							if(!Smits.PhyloCanvas.Render.Style[i]) {  Smits.PhyloCanvas.Render.Style[i] = {}; };
							for(var j in styles[i][0]){
								if(j != '_attributes' && j != '_children' && j != 'type'){
									Smits.PhyloCanvas.Render.Style[i][j.replace('_', '-')] = styles[i][0][j];		// This is quite painful, as xml does not allow dashes
								}
							}
						}
						
					}
				}
			}
			
			// Custom Parameters
			if(render && render.parameters){
				recursiveProcessParameters(render.parameters[0]);
			}			
			
			// Charts
			if(render && render.charts){
				var charts = render.charts[0];
				for (var i in charts){
					if(i != '_children' && i != 'Text'){
						for(var j in charts[i]){
							if(charts[i][j].type == "binary"){
								charts[i][j].chart = i;
								Smits.PhyloCanvas.Render.Parameters.binaryCharts.push(charts[i][j]);
							} else if (charts[i][j].type == "integratedBinary"){
								charts[i][j].chart = i;
								Smits.PhyloCanvas.Render.Parameters.integratedBinaryCharts.push(charts[i][j]);								
							} else if (charts[i][j].type == "bar"){
								charts[i][j].chart = i;
								Smits.PhyloCanvas.Render.Parameters.barCharts.push(charts[i][j]);
							}
						}
					}
				}
			}			
			
		}
		
		root = recursiveProcessRoot(root);

	}

}();

Smits.PhyloCanvas.PhyloxmlParse.prototype = {

};Smits.PhyloCanvas.NexmlParse = function(){

	var mLevel = 0,
	mNewickLen = 0,
	root,
	validate,
	nexEdges,
	nexNodes,
		
	recursiveParse = function(nexnode, nexlen, parentNode){
		var node = new Smits.PhyloCanvas.Node();
		if(parentNode){
			node.level = parentNode.level + 1;
		}
		
		for(var i = 0; i < nexEdges.length; i++){
			if(nexEdges[i].source == nexnode.id){
				for(var j = 0; j < nexNodes.length; j++){
					if(nexEdges[i].target == nexNodes[j].id){
						node.children.push(recursiveParse(nexNodes[j], nexEdges[i].length, node));
					}
				}
			}
		}

		if(node && node.level > 0 && !node.len){
			node.len = 1;
		} 
		
		if(nexlen) {
			node.len = Smits.Common.roundFloat(nexlen, 4);			// round to 4 decimal places
			if(node.len == 0){
				node.len = 0.0001;
			}			
		}

		if(nexnode.label){
			node.type = 'label';
			node.name = nexnode.label;
			if(nexnode.style){
				node.style = nexnode.style;
			}
		} 


			
		return node;
	},
	
	recursiveProcessRoot = function(node, parentNode){
		if(node.children && node.children.length){
			for( var i = 0; i < node.children.length; i++){
				var child = node.children[i];
				child.newickLen = Math.round( (child.len + node.newickLen) *10000)/10000;
				if(child.level > mLevel) mLevel = child.level;
				if(child.newickLen > mNewickLen) mNewickLen = child.newickLen;
				if(child.children.length > 0){
					recursiveProcessRoot(child, node); 
				}				
			}
		}
		return node;
	},
	recursiveProcessParameters = function(parametersEl, treeType){
		for (var i in parametersEl){
			if(i != '_children' && i != 'Text'){
				if(i == 'rectangular' || i == 'circular'){
					recursiveProcessParameters(parametersEl[i][0], i);
				} else {
					if(!Smits.PhyloCanvas.Render.Parameters[i]) {  Smits.PhyloCanvas.Render.Parameters[i] = {}; };
					Smits.PhyloCanvas.Render.Parameters.set(i, parametersEl[i][0].Text, treeType);
				}
			}
		}
		return;
	};

	return function(jsonString, inputFormat){
		/* Privileged Methods */
		this.getRoot = function(){
			return root;
		};
		this.getLevels = function(){
			return mLevel;
		};
		this.getNewickLen = function(){
			return mNewickLen;
		};		
		this.getValidate = function(){
			return validate;
		};
		
		
		if(inputFormat.tree && jsonString.trees[0] && jsonString.trees[0].tree[(inputFormat.tree-1)]){
			nexEdges = jsonString.trees[0].tree[(inputFormat.tree-1)].edge;
			nexNodes = jsonString.trees[0].tree[(inputFormat.tree-1)].node;
		} else {
			nexEdges = jsonString.trees[0].tree[0].edge;
			nexNodes = jsonString.trees[0].tree[0].node;
		}

		// Determine Root
		// If defined, default to that node, else use RAV's implementation.
		// RAV 05-22-2011. 
		// It is more robust to search for the root by the tree topology
		// then by looking for a @root attribute. Valid NeXML tree structures
		// always have one node without normal edges pointing into it. The
		// root attribute is used to indicate that this tree is actually rooted.
		// Compare this with nexus/newick: newick strings are always implicitly
		// rooted, even if the tree is called a utree or the [&U] token is used.
		for(var i = 0; i < nexNodes.length; i++) {
			var targetCount = 0;
			if(nexNodes[i].root && nexNodes[i].root == "true"){
				root = nexNodes[i];
				break;
			}
			for(var j = 0; j < nexEdges.length; j++) {
				if(nexEdges[j].target == nexNodes[i].id) {
					targetCount++;
				}
			}
			if ( targetCount == 0 ) {
				root = nexNodes[i];                                     
				break;
			}
		}

		if(root){
			root = recursiveParse(root);
		
			root = recursiveProcessRoot(root);
		} else {
			validate = 'Error. Currently, only rooted NeXML trees are supported.';
		}

	};

}();

Smits.PhyloCanvas.NexmlParse.prototype = {

};Smits.PhyloCanvas.NexmlJsonParse = function(){

	var mLevel = 0,
	mNewickLen = 0,
	root,
	validate,
	nexEdges = [], nexNodes = [],
		
	recursiveParse = function(nexnode, nexlen, parentNode){
		var node = new Smits.PhyloCanvas.Node();
		if(parentNode){
			node.level = parentNode.level + 1;
		}
		
		for(var i = 0; i < nexEdges.length; i++){
			if(nexEdges[i].source == nexnode.id){
				for(var j = 0; j < nexNodes.length; j++){
					if(nexEdges[i].target == nexNodes[j].id){
						node.children.push(recursiveParse(nexNodes[j], nexEdges[i].length, node));
					}
				}
			}
		}
		
		if(nexlen) {
			node.len = Smits.Common.roundFloat(nexlen, 4);			// round to 4 decimal places
			if(node.len == 0){
				node.len = 0.0001;
			}			
		}

		if(nexnode.label){
			node.type = 'label';
			node.name = nexnode.label;
			if(nexnode.accession){
				node.accession = nexnode.accession;
			}
			if(nexnode.style){
				node.style = nexnode.style;
			}
			if(nexnode.bgStyle){
				node.bgStyle = nexnode.bgStyle;
			}
		} 

		if(nexnode.chart){
			node.chart = nexnode.chart;
		}

		// Validation
		if(node && node.level > 1){
			if(!node.len){
				validate = 'Error. Please include Branch Lengths - we only draw rooted phylogenetic trees.';
			}
		} 
			
		return node;
	},
	
	recursiveProcessRoot = function(node, parentNode){
		if(node.children && node.children.length){
			for( var i = 0; i < node.children.length; i++){
				var child = node.children[i];
				child.newickLen = Math.round( (child.len + node.newickLen) *10000)/10000;
				if(child.level > mLevel) mLevel = child.level;
				if(child.newickLen > mNewickLen) mNewickLen = child.newickLen;
				if(child.children.length > 0){
					recursiveProcessRoot(child, node); 
				}				
			}
		}
		return node;
	},
	
	recursiveProcessParameters = function(parametersEl, treeType){
		for (var i in parametersEl){
			if(i != '_children' && i != 'Text'){
				if(i == 'rectangular' || i == 'circular'){
					recursiveProcessParameters(parametersEl[i], i);
				} else {
					if(!Smits.PhyloCanvas.Render.Parameters[i]) {  Smits.PhyloCanvas.Render.Parameters[i] = {}; };
					Smits.PhyloCanvas.Render.Parameters.set(i, parametersEl[i], treeType);
				}
			}
		}
		return;
	};

	return function(inputFormat){
		/* Privileged Methods */
		this.getRoot = function(){
			return root;
		};
		this.getLevels = function(){
			return mLevel;
		};
		this.getNewickLen = function(){
			return mNewickLen;
		};		
		this.getValidate = function(){
			return validate;
		};
		
		var jsonString = inputFormat.nexmlJson.nexml;
		

		/* RENDER STYLES */
		var render = jsonString.render;
		
		// Custom Styles
		if(render && render.styles){
			var styles = render.styles;
			for (var i in styles){
				if(i != '_children' && i != 'Text'){
					if(styles[i]['@attributes']['type'] && styles[i]['@attributes']['type'] == "radialGradient" && Raphael.svg){
						// radialGradient only supported by SVG
						styles[i]['name'] = i;
						styles[i]['type'] = styles[i]['@attributes']['type'];
						Smits.PhyloCanvas.Render.Style[i] = styles[i];
						if(!Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList']) { Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'] = [] };
						Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'].push(i); 
					} else {
						if(!Smits.PhyloCanvas.Render.Style[i]) {  Smits.PhyloCanvas.Render.Style[i] = {}; };
						for(var j in styles[i]['@attributes']){
							if(j != '_attributes' && j != '_children' && j != 'type'){
								Smits.PhyloCanvas.Render.Style[i][j.replace('_', '-')] = styles[i]['@attributes'][j];		// This is quite painful, as xml does not allow dashes
							}
						}
					}
				}
			}
		}
		// Custom Parameters
		if(render && render.parameters){
			recursiveProcessParameters(render.parameters);
		}
		
		// Charts
		if(render && render.charts){
			var charts = render.charts;
			for (var i in charts){
				charts[i]['@attributes'].chart = i;
				if(charts[i]['@attributes'].type == "binary"){
					Smits.PhyloCanvas.Render.Parameters.binaryCharts.push(charts[i]['@attributes']);
				} else if(charts[i]['@attributes'].type == "integratedBinary"){
					Smits.PhyloCanvas.Render.Parameters.integratedBinaryCharts.push(charts[i]['@attributes']);
				} else if(charts[i]['@attributes'].type == "bar"){
					Smits.PhyloCanvas.Render.Parameters.barCharts.push(charts[i]['@attributes']);
				}
			}
		}				
			
		if(inputFormat.tree && jsonString.trees[0] && jsonString.trees[0].tree[(inputFormat.tree-1)]){
			nexEdges = jsonString.trees[0].tree[(inputFormat.tree-1)].edge;
			nexNodes = jsonString.trees[0].tree[(inputFormat.tree-1)].node;
		} else {
			for(var i = 0; i < jsonString.trees.tree.edge.length; i++){
				nexEdges.push(jsonString.trees.tree.edge[i]['@attributes']);
			}
			for(var i = 0; i < jsonString.trees.tree.node.length; i++){
				var node = jsonString.trees.tree.node[i]['@attributes'];
				if(node.label){
					node.chart = jsonString.trees.tree.node[i].chart;
				}
				nexNodes.push(node);
			}
		}
		
		for(var i = 0; i < nexNodes.length; i++){
			if(nexNodes[i].root && nexNodes[i].root == "true"){
				root = nexNodes[i];
			}
		}
		
		if(root){
			root = recursiveParse(root);
		
			root = recursiveProcessRoot(root);
		} else {
			validate = 'Error. Currently, only rooted NeXML trees are supported.';
		}

	};

}();

Smits.PhyloCanvas.NexmlParse.prototype = {

};Smits.PhyloCanvas.Render = {};Smits.PhyloCanvas.Render.Style = {

	/* Default Styles */
	
	line: {
		"stroke": 		'rgb(0,0,0)',
		"stroke-width":	1
	},
	
	text: {
		"font-family":	'Verdana',
		"font-size":	12,
		"text-anchor":	'start'
	},
	
	path: {
		"stroke": 		'rgb(0,0,0)',
		"stroke-width":	1	
	},
	
	connectedDash : {
		"stroke": 			'rgb(200,200,200)',
		"stroke-dasharray":	". "
	},
	
	textSecantBg : {
		"fill": 	'#EEE',
		"stroke":	'#DDD'
	},
	
	highlightedEdgeCircle : {
		"fill": 	'red'
	},
	
	barChart : {
		fill:		'#003300',
		stroke:		'#DDD'
	},
	
	getStyle : function(requestStyle, fallbackStyle){
		if(this[requestStyle]){
			return this[requestStyle];
		} else {
			return this[fallbackStyle];
		}
	
	}



};Smits.PhyloCanvas.Render.Parameters = {

	/* DEFAULT PARAMETERS */
	jsOverride: 0,				// If set, js will override chart's file setting
	
	/** Phylogram parameters are separated because they behave very differently **/
	
	/* Rectangular Phylogram */
	Rectangular : {
		bufferX			: 200, 			// Reduces the available canvas space for tree branches, allowing
										// for more space for the textual/charting components
		paddingX		: 10,
		paddingY		: 20,
		bufferInnerLabels : 10, 		// Pixels
		bufferOuterLabels : 5, 			// Pixels
		minHeightBetweenLeaves : 10,  	// Should probably set pretty low, as clipping may occur if it needs to be implemented		
		
		alignPadding	: 0,			// Pixels to push the labels out by - this extension should be 
										// compensated by an increase in bufferX too
		alignRight		: false,
		
		showScaleBar	: false			// (STRING,  e.g. "0.05") Shows a scale bar at the bottom of the tree
	},
	
	/* Circular Phylogram */
	Circular : {
		bufferRadius 		: 0.33,		// Margins of Tree Circle
										// If > 1, it is in pixels
										// If < 1, it is a percentage of the full canvas size		
		bufferAngle 		: 20,		// controls split size in circle		
		initStartAngle 		: 160,		
		innerCircleRadius 	: 0,
		minHeightBetweenLeaves : 5,

		/* Labels */
		bufferInnerLabels : 2, 			// Pixels
		bufferOuterLabels : 5 			// Pixels
	},
	
	/* Charts */
	binaryCharts : [],
	integratedBinaryCharts : [],
	barCharts : [],

		/* Binary Defaults */
		binaryChartBufferInner : 5, 
		binaryChartBufferSiblings : 0.01,
		binaryChartThickness : 15,
		binaryChartDisjointed : false,
			
		/* Bar Defaults */
		barChartBufferInner : 3,
		barChartHeight : 50,
		barChartWidth : 0.5,	// If > 1, it is in pixels
								// If < 1, it is a percentage of the node width 
						
		/* 
			Rollover Events 
				At minimum, the params object has the following properties:
					.svg
					.node
					.x
					.y
					.textEl
		*/
		mouseRollOver : function(params) {
			if(params.node.edgeCircleHighlight){
				params.node.edgeCircleHighlight.show();
			} else {
				var circleObject = params.svg.draw(
					new Smits.PhyloCanvas.Render.Circle(
						params.x, params.y, 5,
						{ attr: Smits.PhyloCanvas.Render.Style.highlightedEdgeCircle }
					)
				);
				params.node.edgeCircleHighlight = circleObject[0];
			}					
			params.textEl.attr({ fill: 'red' });
		},
		mouseRollOut : function(params) {
			params.node.edgeCircleHighlight.hide();
			params.textEl.attr({ fill: '#000' });
		},

	set : function(param, value, treeType){
		if(!this.jsOverride){
			if(treeType){
				if(treeType == 'circular'){				
					this['Circular'][param] = parseFloat(value);
				} else if (treeType == 'rectangular'){
					this['Rectangular'][param] = parseFloat(value);
				}
			} else {
				this[param] = parseFloat(value);
			}
		}
	}
};Smits.PhyloCanvas.Render.Line = function(){

	return function(x1, x2, y1, y2, params){
		/* Defaults */	
		this.type = 'line';
		this.attr = Smits.PhyloCanvas.Render.Style.line;
		
		this.x1 = x1;
		this.x2 = x2;
		this.y1 = y1;
		this.y2 = y2;
		
		if(params) {
			Smits.Common.apply(this, params);
			if(params.attr) this.attr = params.attr;
		}

	}
}();Smits.PhyloCanvas.Render.Text = function(){

	return function(x, y, text, params){
		/* Defaults */
		this.type = 'text';
		this.attr = Smits.PhyloCanvas.Render.Style.text;
		
		this.x = x;
		this.y = y;
		this.text = text;
		
		if(params) {
			Smits.Common.apply(this, params);
			if(params.attr) this.attr = params.attr;
		}
	}
}();Smits.PhyloCanvas.Render.Path = function(){
	var attr = Smits.PhyloCanvas.Render.Style.path;
	
	return function(path, params){
		/* Defaults */
		this.type = 'path';
		this.attr = Smits.PhyloCanvas.Render.Style.path;
		
		this.path = path;
		if(params) {
			Smits.Common.apply(this, params);
			if(params.attr) this.attr = params.attr;
		}

	}
}();Smits.PhyloCanvas.Render.Circle = function(){

	return function(x, y, radius, params){
		/* Defaults */	
		this.type = 'circle';
	
		this.x = x;
		this.y = y;
		this.radius = radius;
		
		if(params) {
			Smits.Common.apply(this, params);
			if(params.attr) this.attr = params.attr;
		}
		
	}
}();Smits.PhyloCanvas.Render.SVG = function(){
	var divId,
		canvasSize;
		
	return function(sDivId, canvasWidth, canvasHeight){
	
		/* CONSTRUCTOR */
		divId = sDivId;
		this.canvasSize = [canvasWidth, canvasHeight];
		
		this.svg = Raphael(sDivId, this.canvasSize[0], this.canvasSize[1]);
		
	}
	
}();

Smits.PhyloCanvas.Render.SVG.prototype = {

	render : function(){
		var instructs = this.phylogramObject.getDrawInstructs();
		console.log('render', this.phylogramObject.getDrawInstructs());
		for (var i = 0; i < instructs.length; i++) {
		   if(instructs[i].type == 'line'){
				var line = this.svg.path(["M", instructs[i].x1, instructs[i].y1, "L", instructs[i].x2, instructs[i].y2]).attr(Smits.PhyloCanvas.Render.Style.line);
			} else if(instructs[i].type == 'path'){
				var path = this.svg.path(instructs[i].path).attr(instructs[i].attr);			
			} else if(instructs[i].type == 'circle'){
				var path = this.svg.circle(instructs[i].x, instructs[i].y, instructs[i].radius).attr({
					"stroke": 'red'
				});
			} else {
				var text = this.svg.text(instructs[i].x, instructs[i].y, instructs[i].text).attr(Smits.PhyloCanvas.Render.Style.text);
				if(instructs[i].attr){
					text.attr(instructs[i].attr);
				}
				if(instructs[i].rotate){
					text.rotate(instructs[i].rotate);
				}
				
				var bbox = text.getBBox();
				var hyp = Math.sqrt( (bbox.height * bbox.height) + (bbox.width * bbox.width) );	// get hypotenuse
				
			} 
		}
	},
	
	draw : function(instruct){
		var obj, 
			param;

	   if(instruct.type == 'line'){
			obj = this.svg.path(["M", instruct.x1, instruct.y1, "L", instruct.x2, instruct.y2]).attr(Smits.PhyloCanvas.Render.Style.line);
		} else if(instruct.type == 'path'){
			obj = this.svg.path(instruct.path).attr(instruct.attr);			
		} else if(instruct.type == 'circle'){
			obj = this.svg.circle(instruct.x, instruct.y, instruct.radius).attr({
				"stroke": 'red'
			});
		} else if(instruct.type == 'text'){
			obj = this.svg.text(instruct.x, instruct.y, instruct.text).attr(Smits.PhyloCanvas.Render.Style.text);
			if(instruct.attr){
				obj.attr(instruct.attr);
			}
			if(instruct.rotate){
				obj.rotate(instruct.rotate);
			}
			
			var bbox = obj.getBBox();
			param = Math.sqrt( (bbox.height * bbox.height) + (bbox.width * bbox.width) );	// get hypotenuse
		} 

		return [obj, param];
	}

};Smits.PhyloCanvas.Render.Phylogram = function(){

	var svg,
	sParams = Smits.PhyloCanvas.Render.Parameters.Rectangular, 	// Easy Reference
	canvasX, canvasY,
	scaleX, scaleY, maxBranch,
	minHeightBetweenLeaves,
	firstBranch = true,
	absoluteY = 0, maxLabelLength = 0,
	outerX, outerY, outerRadius,
	x1, x2, y1, y2, 
	positionX, positionY,
	bufferX, paddingX, paddingY, labelsHold = [],
	
	textPadding = function (y){
		return y + Math.round(y / 4);
	},
	
	rectLinePathArray = function (x1, y1, x2, y2){
		return ["M", x1, y1, "L", x2, y1, "L", x2, y2, "L", x1, y2, "Z"];
	},
	
	recursiveCalculateNodePositions = function (node, positionX){
		if(node.len && firstBranch == false && node.children.length == 0){ 
			absoluteY = Smits.Common.roundFloat(absoluteY + scaleY, 4);
		}
		
		if(node.children.length > 0){
			var nodeCoords = [], x1,x2,y1,y2;
			if(node.len){ // draw stem
				x1 = positionX;
				x2 = positionX = Smits.Common.roundFloat(positionX + (scaleX * node.len), 4);
				y1 = absoluteY + (node.getMidbranchPosition(firstBranch) * scaleY);
				y2 = y1;
				svg.draw(new Smits.PhyloCanvas.Render.Line(x1, x2, y1, y2));
			}
			
			if(node.name){ // draw bootstrap values
				var attr = {};
				attr = Smits.PhyloCanvas.Render.Style.getStyle('bootstrap', 'text');
				if(node.uri) { attr.href = node.uri };
				if(node.description) {attr.title = node.description };
				if(node.level == 0){ 
					var innerY2 = absoluteY + (node.getMidbranchPosition(firstBranch) * scaleY);
				} else {
					var innerY2 = y2;
				}
				
				svg.draw(
					new Smits.PhyloCanvas.Render.Text(
						(x2 || positionX) + 5, innerY2,
						node.name,
						{
							attr: attr
						}
					)
				);			
			}
			
			if(node.children && node.children.length){
				for(var i = 0; i < node.children.length; i++){
					var child = node.children[i];
					nodeCoords.push(recursiveCalculateNodePositions(child, positionX));
				}
			}
			
			var flatNodeCoords = []; // establish vertical bounds
			for ( var i = 0; i < nodeCoords.length; i++ ){
				if(nodeCoords[i][0]) flatNodeCoords.push(nodeCoords[i][0]);
				if(nodeCoords[i][1]) flatNodeCoords.push(nodeCoords[i][1]);
			}
			var verticalY1 = Math.min.apply(null, flatNodeCoords );
			var verticalY2 = Math.max.apply(null, flatNodeCoords);
			
			// draw vertical
			// hack: little elbows at ends in order to prevent stair-effects at edges
			svg.draw( 
				new Smits.PhyloCanvas.Render.Path( 
					[
						"M", positionX + 0.0001, verticalY1,
						"L", positionX, verticalY1,
						"L", positionX, verticalY2,
						"L", positionX + 0.0001, verticalY2
					],
					{ attr : Smits.PhyloCanvas.Render.Style.line }
				)				
			);
			
		} else {
			// label
			x1 = positionX;
			x2 = Smits.Common.roundFloat(positionX + (scaleX * node.len), 2);
			y1 = absoluteY;
			y2 = absoluteY;
				
			// preserve for later processing
			node.y = absoluteY;
			labelsHold.push(node);				
				
			svg.draw(new Smits.PhyloCanvas.Render.Line(x1, x2, y1, y2));
			if(sParams.alignRight){
				svg.draw(
					new Smits.PhyloCanvas.Render.Path(
						["M", x2, y1, "L", sParams.alignPadding + maxBranch, y2],
						{ attr : Smits.PhyloCanvas.Render.Style.connectedDash }
					)
				);			
			}
			
			if(node.name){
				var attr = {};
				if(node.style){
					attr = Smits.PhyloCanvas.Render.Style.getStyle(node.style, 'text');
				}
				attr["text-anchor"] = 'start';
				if(node.uri) { attr.href = node.uri };
				if(node.description) {attr.title = node.description };
				
				var draw = svg.draw(
					new Smits.PhyloCanvas.Render.Text(
						sParams.alignRight ? maxBranch + sParams.bufferInnerLabels + sParams.alignPadding : x2 + sParams.bufferInnerLabels, y2,
						node.name,
						{
							attr: attr
						}
					)
				);				
				maxLabelLength = Math.max(draw[1], maxLabelLength);
				

				// Rollover, Rollout and Click Events
				if(Smits.PhyloCanvas.Render.Parameters.mouseRollOver){
					Smits.Common.addRaphEventHandler(
						draw[0], 
						'mouseover', 
						Smits.PhyloCanvas.Render.Parameters.mouseRollOver, 
						{ svg: svg, node: node, x: x2, y: y2, textEl: draw[0] }
					);
				}
				if(Smits.PhyloCanvas.Render.Parameters.mouseRollOut){
					Smits.Common.addRaphEventHandler(
						draw[0], 
						'mouseout', 
						Smits.PhyloCanvas.Render.Parameters.mouseRollOut, 
						{ svg: svg, node: node, x: x2, y: y2, textEl: draw[0] }
					);				
				}
				if(Smits.PhyloCanvas.Render.Parameters.onClickAction){
					Smits.Common.addRaphEventHandler(
						draw[0], 
						'click', 
						Smits.PhyloCanvas.Render.Parameters.onClickAction, 
						{ svg: svg, node: node, x: x2, y: y2, textEl: draw[0] }
					);				
				}
			}
			
			if(firstBranch){
				firstBranch = false;
			}
		
		}
		
		return [y1, y2];

	},
	
	drawScaleBar = function (){
		y = absoluteY + scaleY;
		x1 = 0;
		x2 = sParams.showScaleBar * scaleX;
		svg.draw(new Smits.PhyloCanvas.Render.Line(x1, x2, y, y));
		svg.draw(new Smits.PhyloCanvas.Render.Text(
			(x1+x2)/2, 
			y-8, 
			sParams.showScaleBar)
		);
	},
	
	renderBinaryChart = function(x, groupName, params){
		var bufferInner = (params && params.bufferInner ? params.bufferInner : 0) | Smits.PhyloCanvas.Render.Parameters.binaryChartBufferInner,
			bufferSiblings = (params && params.bufferSiblings ? params.bufferSiblings * scaleY : 0) | (Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings < 1 ? scaleY * Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings : Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings),		
			thickness = (params && params.thickness ? params.thickness : 0) | Smits.PhyloCanvas.Render.Parameters.binaryChartThickness,
			beginY;
			
		for(var i = 0; i < labelsHold.length; i++){
			var node = labelsHold[i];
			svg.draw(
				new Smits.PhyloCanvas.Render.Path(
					rectLinePathArray(
						x + bufferInner,
						node.y - (scaleY/2) + (bufferSiblings/2),
						x + bufferInner + thickness, 
						node.y + (scaleY/2) - (bufferSiblings/2)
					),
					{ attr: Smits.PhyloCanvas.Render.Style.getStyle(node.chart[groupName], 'textSecantBg') }
				)
			);			
		}
		return x + bufferInner + thickness;
	},
	
	renderBarChart = function(x, groupName, params){
		var allValues = [], maxValue,
			bufferInner = params && params.bufferInner ? params.bufferInner : 0 | Smits.PhyloCanvas.Render.Parameters.barChartBufferInner,
			height = params && params.height ? params.height : 0 | Smits.PhyloCanvas.Render.Parameters.barChartHeight,
			width = params && params.width ? (params.width < 1 ? scaleY * params.width : params.width ) : 0 | (Smits.PhyloCanvas.Render.Parameters.barChartWidth < 1 ? scaleY * Smits.PhyloCanvas.Render.Parameters.barChartWidth : Smits.PhyloCanvas.Render.Parameters.barChartWidth),
			scaleHeight = 0;
		
		// Need to get max value
		for(var i = 0; i < labelsHold.length; i++){
			allValues.push(labelsHold[i].chart[groupName]);
		}
		maxValue = Math.max.apply(null, allValues);
		scaleHeight = Smits.Common.roundFloat(height / maxValue, 4);
		
		for(var i = 0; i < labelsHold.length; i++){
			var node = labelsHold[i];
			svg.draw(
					new Smits.PhyloCanvas.Render.Path(
						rectLinePathArray(
							x + bufferInner,
							node.y - (width/2),
							x + bufferInner + (scaleHeight * node.chart[groupName]), 
							node.y + (width/2)
						),					
						{ attr: Smits.PhyloCanvas.Render.Style.getStyle(node.chart[groupName], 'barChart') }
					)
				);					
		}
		
		return x + bufferInner + height;
	};
	
	return function(sSvg, dataObject){	
		/* Privileged Methods */
		this.getCanvasSize = function(){
			return [canvasX, canvasY];
		};		
		this.getRoot = function(){
			return dataObject.getRoot();
		};
		
		/* CONSTRUCTOR */
		if(dataObject.getValidate()){   // Validate
			svg.draw(0,0, dataObject.getValidate());
		}
		
		svg = sSvg;
		var node = dataObject.getRoot();
		var mNewickLen = dataObject.getNewickLen();
		
		canvasX = svg.canvasSize[0];			// Full Canvas Width
		canvasY = svg.canvasSize[1];			// Full Canvas Height
		
		bufferX = sParams.bufferX;
		paddingX = sParams.paddingX;
		paddingY = sParams.paddingY;
		minHeightBetweenLeaves = sParams.minHeightBetweenLeaves;

		absoluteY = paddingY;
		
		scaleX = Math.round((canvasX - bufferX - paddingX*2) / mNewickLen);
		scaleY = Math.round((canvasY - paddingY*2) / (sParams.showScaleBar ? node.getCountAllChildren() : node.getCountAllChildren() - 1 ) );
		if(scaleY < minHeightBetweenLeaves){
			scaleY = minHeightBetweenLeaves;
		}
		maxBranch = Math.round( canvasX - bufferX - paddingX*2 );	
		
		if(Smits.PhyloCanvas.Render.Parameters.binaryCharts.length || Smits.PhyloCanvas.Render.Parameters.barCharts.length){
			sParams.alignRight = true;
		}
		
		recursiveCalculateNodePositions(node, paddingX);
		
		// Draw Scale Bar
		if(sParams.showScaleBar){
			drawScaleBar();
		}
		
		outerX = maxBranch + maxLabelLength + sParams.bufferInnerLabels;
		// Draw secant highlights
		if(Smits.PhyloCanvas.Render.Parameters.binaryCharts.length){
			var binaryCharts = Smits.PhyloCanvas.Render.Parameters.binaryCharts;
			for(var i in binaryCharts){
				outerX = renderBinaryChart(outerX, binaryCharts[i].chart, binaryCharts[i]);
			}
		}		
		
		// Draw Bar Chart
		if(Smits.PhyloCanvas.Render.Parameters.barCharts.length){
			var barCharts = Smits.PhyloCanvas.Render.Parameters.barCharts;
			for(var i in barCharts){
				outerRadius = renderBarChart(outerX, barCharts[i].chart, barCharts[i]);
			}
		}				

	}
}();

Smits.PhyloCanvas.Render.Phylogram.prototype = {

};Smits.PhyloCanvas.Render.CircularPhylogram = (function(){

	var svg,
		sParams = Smits.PhyloCanvas.Render.Parameters.Circular, 	// Easy Reference
		canvasX, canvasY, canvasMinEdge,
		scaleRadius, scaleAngle,
		minHeightBetweenLeaves,
		innerCircleRadius,
		firstBranch = true,
		absoluteY = 0, cx, cy, maxBranch, 
		labelsHold = [], bgLabelsHold = [],
		bufferRadius, bufferAngle, outerRadius,
		maxLabelLength = 0,
		initStartAngle,
		rad = (Math.PI / 180);

	function secPosition(r, deg){
		deg += initStartAngle;
		return [ 
			Smits.Common.roundFloat(cx + r * Math.sin(deg * rad), 4), 
			Smits.Common.roundFloat(cy + r * Math.cos(deg * rad), 4)
		]; // x,y
	};
	function rotateTextByY(yCoord){
		var rotateAngle = normalizeAngle( 90 - yCoord - initStartAngle );
			
		if(rotateAngle > 90 && rotateAngle < 270){
			rotateAngle += 180;
			var alignment = "end";
		} else {
			var alignment = "start";
		}	
		
		return [rotateAngle, alignment];
	};
	function secant(r, startAngle, endAngle, params){
		var startPos = secPosition(r, startAngle);
		var endPos = secPosition(r, endAngle);
		var arr = [],
			n, inv = 0;
		
		if(Math.abs(normalizeAngle(endAngle-startAngle)) > 180) {
			n = 1;
		} else {
			n = -1;
		}
		
		// Parameter changes
		if(params && params.invertSecant){
			n *= -1;
			inv = 1;
		}
		if(params && params.noMove){
		} else {
			arr.push('M');
		}
		
		arr.push(startPos[0], startPos[1], "A", r, r, 0, n < 1 ? 0 : 1, inv, endPos[0], endPos[1]);
		return arr;
	};
	function secLinePath(deg, x1, x2, params){
		var arr = [];
		var startPos = secPosition(x1, deg);
		var endPos = secPosition(x2, deg);
		if(params && params.noMove){
		} else {
			arr.push('M');
		}
		arr.push(startPos[0], startPos[1], "L", endPos[0], endPos[1]);
		return arr;
	};
	function normalizeAngle(ang){
		while(ang > 360 || ang < 0){
			if(ang > 360){
				ang -= 360;
			} else if (ang < 0){
				ang += 360;
			}
		}
		return ang;
	};
	function sector(r1, r2, y1, y2){
		if(!r2 && r1.length > 1){
			var y2 = r1[3];
			var y1 = r1[2];
			var r2 = r1[1];
			var r1 = r1[0];
		}
		var arr = array_merge( "M",
			secant(
				r1, 
				y1,
				y2, 
				{ noMove: 1, invertSecant: 0}
			), "L",
			secant(
				r2, 
				y2, 
				y1, 
				{ noMove: 1, invertSecant: 1}
			),
			'Z'
		);	
		return arr;
	};
	
	function recursiveCalculateNodePositions(node, positionX){
		positionX = positionX;

		if(node.len){ // If first branch, pad only margin
			if(firstBranch){
				absoluteY = bufferAngle || 1;		// Has to be at least 1
				
			} else {
				if(node.children.length == 0) absoluteY = Smits.Common.roundFloat(absoluteY + scaleAngle, 4);
			}
		}
		if(node.children.length > 0){
			var nodeCoords = [], x1,x2,y1,y2;
			x1 = positionX;
			x2 = positionX += Smits.Common.roundFloat(scaleRadius * node.len, 4);
				
		
			if(node.name){ // draw bootstrap values
				
			}
			
			if(node.children && node.children.length){
				for(var i = 0; i < node.children.length; i++){
					var child = node.children[i];
					var y = recursiveCalculateNodePositions(child, positionX);
					if(y > 0) nodeCoords.push(y);					
				}
			}
			
			var minAngle = Smits.Common.roundFloat(Math.min.apply(null, nodeCoords ), 4);
			var maxAngle = Smits.Common.roundFloat(Math.max.apply(null, nodeCoords ), 4);
			
			// hack: little elbows at ends in order to prevent stair-effects at edges
			if(node.level != 0){
				svg.draw(
					new Smits.PhyloCanvas.Render.Path(
						array_merge(
							"M", secPosition(positionX + 0.01, minAngle), 
							"L", secant(positionX, minAngle, maxAngle, {noMove: true}),
							"L", secPosition(positionX + 0.01, maxAngle)
						)
					)
				);
			}
			
			if(node.len){ // draw stem
				y1 = Smits.Common.roundFloat( minAngle + (maxAngle-minAngle)/2, 4 );
				svg.draw(new Smits.PhyloCanvas.Render.Path(secLinePath(y1, x1, x2)));
			}			
			
		} else {			
			// LABEL
			
			// preserve for later processing
			node.y = absoluteY;
			labelsHold.push(node);
			
			x1 = positionX;
			x2 = positionX =  Smits.Common.roundFloat(positionX + (scaleRadius * node.len));
			y1 = absoluteY;
				
			svg.draw(new Smits.PhyloCanvas.Render.Path(secLinePath(y1, x1, x2)));
			svg.draw(
				new Smits.PhyloCanvas.Render.Path(
					secLinePath(y1, x2, maxBranch), 
					{ attr : Smits.PhyloCanvas.Render.Style.connectedDash }
				)
			);
			
			
			
			if(node.name){
				var pos = secPosition(maxBranch + sParams.bufferInnerLabels, y1);
				var rotateParam = rotateTextByY(y1);
				var rotateAngle = rotateParam[0];
				var alignment = rotateParam[1];

				var attr = {};
				if(node.style){
					Smits.Common.apply(attr, Smits.PhyloCanvas.Render.Style.getStyle(node.style, 'text'));
				}
				attr["text-anchor"] = alignment;
				if(node.uri) { attr.href = node.uri };
				if(node.description) {attr.title = node.description };
				
				var draw = svg.draw(
					new Smits.PhyloCanvas.Render.Text(
						pos[0], pos[1], 
						node.name,
						{
							attr: attr,
							rotate: [rotateAngle, pos[0], pos[1]]
						}
					)
				);
				
				// Background Style
				if(node.bgStyle){
					bgLabelsHold.push([node.bgStyle, y1]);
				}
				
				// Rollover, Rollout and Click Events
				var pos = secPosition(x2, y1);
				if(Smits.PhyloCanvas.Render.Parameters.mouseRollOver){
					Smits.Common.addRaphEventHandler(
						draw[0], 
						'mouseover', 
						Smits.PhyloCanvas.Render.Parameters.mouseRollOver, 
						{ svg: svg, node: node, x: pos[0], y: pos[1], textEl: draw[0] }
					);
				}
				if(Smits.PhyloCanvas.Render.Parameters.mouseRollOut){
					Smits.Common.addRaphEventHandler(
						draw[0], 
						'mouseout', 
						Smits.PhyloCanvas.Render.Parameters.mouseRollOut, 
						{ svg: svg, node: node, x: pos[0], y: pos[1], textEl: draw[0] }
					);				
				}
				if(Smits.PhyloCanvas.Render.Parameters.onClickAction){
					Smits.Common.addRaphEventHandler(
						draw[0], 
						'click', 
						Smits.PhyloCanvas.Render.Parameters.onClickAction, 
						{ svg: svg, node: node, x: pos[0], y: pos[1], textEl: draw[0] }
					);							
				}
				
				maxLabelLength = Math.max(draw[1], maxLabelLength);
			}
		}
		if(firstBranch){
			firstBranch = false;
		}
		return y1;
	};

	
	function array_merge(arr) {
		var merged = arr;
		for (var i = 1; i < arguments.length; i++) {
			merged = merged.concat(arguments[i]);
		}
		return merged;
	};
	
	function renderBackground(){
		var arr = [];
		
		// Highlighted Labels
		if(bgLabelsHold.length > 0){
		
			// Setup Gradients if defined
			if(Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList']){
				for(var i = 0; i < Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'].length; i++){
					var gradientName = Smits.PhyloCanvas.Render.Style['jsphylosvgGradientList'][i];
					var radialEl = Smits.Common.createGradientEl(gradientName, Smits.PhyloCanvas.Render.Style[gradientName], [cx, cy, maxBranch + maxLabelLength + sParams.bufferOuterLabels]);
					svg.svg.defs.appendChild(radialEl);
				}
			}		
			
			for(var i = 0; i < bgLabelsHold.length; i++){
				if(i != bgLabelsHold.length - 1 && bgLabelsHold[i][0] == bgLabelsHold[(i+1)][0]){
					bgLabelsHold[(i+1)][2] = bgLabelsHold[i][2] ? bgLabelsHold[i][2] : bgLabelsHold[i][1];
					continue;
				}
				
				var arr = sector(
					maxBranch, 
					maxBranch + maxLabelLength + sParams.bufferOuterLabels, 
					bgLabelsHold[i][2] ? bgLabelsHold[i][2] - scaleAngle/2 : bgLabelsHold[i][1] - scaleAngle/2, 
					bgLabelsHold[i][1] + scaleAngle/2
				);			
				var attr = Smits.PhyloCanvas.Render.Style.getStyle(bgLabelsHold[i][0], 'textSecantBg');
				var bgObj = svg.draw(
					new Smits.PhyloCanvas.Render.Path(
						arr, 
						{ attr: attr.type ? {} : attr}
					)
				);
				//if(attr.type && attr.type == "radialGradient") { bgObj[0].node.setAttribute('class', 'jsphylosvg-' + attr.name); };
				if(attr.type && attr.type == "radialGradient") { bgObj[0].node.setAttribute('fill', 'url(#' + attr.name + ')'); };
				if(attr.type && attr.type == "radialGradient") { bgObj[0].node.setAttribute('stroke', 'none'); };
				bgObj[0].toBack(); 		// Put it behind the labels
			}
		}
		
		// Neutral Background
		var arr = sector(
			maxBranch, 
			maxBranch + maxLabelLength + sParams.bufferOuterLabels, 
			(bufferAngle || 1) - (scaleAngle/2), 
			360  - (scaleAngle/2)
		);
		var bgObj = svg.draw(
			new Smits.PhyloCanvas.Render.Path(
				arr, 
				{ attr: Smits.PhyloCanvas.Render.Style.textSecantBg }
			)
		);
		
		bgObj[0].toBack(); 		// Put it behind the labels
		
		return maxBranch + maxLabelLength + sParams.bufferOuterLabels;
	};
	
	function renderBinaryChart(outerRadius, groupName, params){
		var bufferInner = (params && params.bufferInner) ? parseFloat(params.bufferInner) : Smits.PhyloCanvas.Render.Parameters.binaryChartBufferInner,
			bufferSiblings = (params && params.bufferSiblings ? params.bufferSiblings * scaleAngle : 0) | (Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings < 1 ? scaleAngle * Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings : Smits.PhyloCanvas.Render.Parameters.binaryChartBufferSiblings),
			thickness = (params && params.thickness) ? parseFloat(params.thickness) : Smits.PhyloCanvas.Render.Parameters.binaryChartThickness,
			disjointed = (params && params.disjointed ? params.disjointed : false) | Smits.PhyloCanvas.Render.Parameters.binaryChartDisjointed,
			isInternal = (params && params.isInternal) ? params.isInternal : false, 
			isFirst = true,
			beginY;
		
		
		for(var i = 0; i < labelsHold.length; i++){
			var node = labelsHold[i];
			if( (!labelsHold[i+1] || node.chart[groupName] !== labelsHold[i+1].chart[groupName] || disjointed) && node.chart[groupName] != "none" ){
				var attr = Smits.PhyloCanvas.Render.Style.getStyle(node.chart[groupName], 'textSecantBg');
				if(isInternal){
					var sectorCoords = [
						maxBranch - bufferInner - thickness,  
						maxBranch - bufferInner,  
						(beginY ? beginY : node.y) - (scaleAngle/2) + (isFirst && !disjointed ? 0 : (bufferSiblings/2)), 
						node.y + (scaleAngle/2) - (i == labelsHold.length-1 && !disjointed ? 0 : (bufferSiblings/2))						
					];				
				} else {
					var sectorCoords = [
						outerRadius + bufferInner,  
						outerRadius + bufferInner + thickness, 
						(beginY ? beginY : node.y) - (scaleAngle/2) + (isFirst && !disjointed ? 0 : (bufferSiblings/2)), 
						node.y + (scaleAngle/2) - (i == labelsHold.length-1 && !disjointed ? 0 : (bufferSiblings/2))
					];
				}
				
				if(attr.label){
					var textAttr = Smits.PhyloCanvas.Render.Style.getStyle(attr.labelStyle, 'text');
					var pos = secPosition( (sectorCoords[0] + sectorCoords[1]) / 2, (sectorCoords[2] + sectorCoords[3]) / 2 );
					var rotateParam = rotateTextByY((sectorCoords[2] + sectorCoords[3]) / 2);
					var rotateLabelBy = normalizeAngle(rotateParam[0] + (textAttr["rotate"] ? parseFloat(textAttr["rotate"]) : 0));

					var rotateAngle = normalizeAngle( 90 - (sectorCoords[2] + sectorCoords[3])/2 - initStartAngle );
					if(rotateAngle > 90 && rotateAngle < 270){
						rotateLabelBy += 180;
					}
					
					if(!textAttr["text-anchor"]){
						textAttr["text-anchor"] = "middle";
					}

					var binText = svg.draw(
						new Smits.PhyloCanvas.Render.Text(
							pos[0], 
							pos[1], 
							attr.label,
							{
								attr: textAttr,
								rotate: rotateLabelBy
							}
						)
					);
					binText[0].toBack();
				}

				if(attr.borderStyle){
					var borderAttr = Smits.PhyloCanvas.Render.Style.getStyle(attr.borderStyle, 'textSecantBg');
					var borderSectorCoords = [
						maxBranch,
						borderAttr.fullsize ? sectorCoords[1] : sectorCoords[0],
						sectorCoords[2],
						sectorCoords[3]
					];
					var binBorder = svg.draw(
						new Smits.PhyloCanvas.Render.Path(
							sector( 
								borderSectorCoords
							),
							{ attr: borderAttr }
						)
					);		
					binBorder[0].toBack();
				}
				
				var binObj = svg.draw(
					new Smits.PhyloCanvas.Render.Path(
						sector( 
							sectorCoords
						),
						{ attr: attr }
					)
				);	
				binObj[0].toBack();
				
				beginY = 0;
				isFirst = false;
			} else {
				if(!beginY){ beginY = node.y; }
				if(node.chart[groupName] == "none"){
					beginY = 0;
				}
			}
			isFirst = false;
		}
		return isInternal ? outerRadius : outerRadius + bufferInner + thickness;
	};
	
	function renderBarChart(outerRadius, groupName, params){
		var allValues = [], maxValue,
			bufferInner = params && params.bufferInner ? parseFloat(params.bufferInner) : Smits.PhyloCanvas.Render.Parameters.barChartBufferInner,
			height = params && params.height ? parseFloat(params.height) : (Smits.PhyloCanvas.Render.Parameters.barChartHeight ? Smits.PhyloCanvas.Render.Parameters.barChartHeight : 0),
			width = params && params.width ? (parseFloat(params.width) < 1 ? scaleAngle * parseFloat(params.width) : parseFloat(params.width) ) : 0 | (Smits.PhyloCanvas.Render.Parameters.barChartWidth < 1 ? scaleAngle * Smits.PhyloCanvas.Render.Parameters.barChartWidth : Smits.PhyloCanvas.Render.Parameters.barChartWidth),
			scaleHeight = 0;
		
		// Need to get max value
		for(var i = 0; i < labelsHold.length; i++){
			allValues.push(labelsHold[i].chart[groupName]);
		}
		maxValue = Math.max.apply(null, allValues);
		scaleHeight = Smits.Common.roundFloat(height / maxValue, 4);
		
		for(var i = 0; i < labelsHold.length; i++){
			var node = labelsHold[i];
			if(node.chart[groupName] > 0){
				svg.draw(
						new Smits.PhyloCanvas.Render.Path(
							sector( 
								outerRadius + bufferInner,  
								outerRadius + bufferInner + (scaleHeight * node.chart[groupName]), 
								node.y - (width/2), 
								node.y + (width/2)
							),
							{ attr: Smits.PhyloCanvas.Render.Style.getStyle(node.chart[groupName], 'barChart') }
						)
				);					
			}
		}
		
		return outerRadius + bufferInner + height;
	};
	
	return function(sSvg, dataObject, bufferRadius){
		/* Privileged Methods */
		this.getCanvasSize = function(){
			return [canvasX, canvasY];
		};
		this.getRoot = function(){
			return dataObject.getRoot();
		};
	
		/* CONSTRUCTOR */
		// Validation
		if(dataObject.getValidate()){   
			sSvg.draw({type: 'text', x: 0, y: sSvg.canvasSize[1] / 3, text: dataObject.getValidate() });
			return
		}				
		
		// Properties Setup
		svg 			= sSvg;
		var node 		= dataObject.getRoot();
		var mNewickLen 	= dataObject.getNewickLen();
		canvasX 		= svg.canvasSize[0];															// Full Canvas Width
		canvasY 		= svg.canvasSize[1];															// Full Canvas Height
		cx 				= canvasX / 2;																	// Set Center Position
		cy 				= canvasY / 2;
		canvasMinEdge 	= Math.min.apply(null, [canvasX,canvasY]);
		
		bufferRadius		= (sParams.bufferRadius > 1) ? sParams.bufferRadius : Smits.Common.roundFloat(canvasMinEdge * sParams.bufferRadius, 4);
		bufferAngle 		= sParams.bufferAngle;							// controls split size in circle		
		innerCircleRadius	= sParams.innerCircleRadius;
		minHeightBetweenLeaves	= sParams.minHeightBetweenLeaves;
		initStartAngle		= sParams.initStartAngle;						// Angle at which the entire tree is rotated
		
		maxBranch			= Math.round( (canvasMinEdge - bufferRadius - innerCircleRadius) / 2);		// maximum branch length
		scaleRadius			= (maxBranch - innerCircleRadius) / mNewickLen;								// scale multiplier to use
		scaleAngle 			= Smits.Common.roundFloat( (360 - bufferAngle) / node.getCountAllChildren(), 4 );		

		// Draw Nodes and Labels
		recursiveCalculateNodePositions(node, innerCircleRadius);
		outerRadius = maxBranch + maxLabelLength + sParams.bufferOuterLabels;

		// Draw integrated secant highlights
		if(Smits.PhyloCanvas.Render.Parameters.integratedBinaryCharts.length){
			var integratedBinaryCharts = Smits.PhyloCanvas.Render.Parameters.integratedBinaryCharts;
			for(var i in integratedBinaryCharts){
				var bufferInner = (integratedBinaryCharts[i].bufferInner ? integratedBinaryCharts[i].bufferInner : Smits.PhyloCanvas.Render.Parameters.binaryChartBufferInner);
				var thickness = (integratedBinaryCharts[i].thickness ? integratedBinaryCharts[i].thickness : Smits.PhyloCanvas.Render.Parameters.binaryChartThickness);
				outerRadius = renderBinaryChart(
					outerRadius - thickness - bufferInner, 
					integratedBinaryCharts[i].chart, 
					integratedBinaryCharts[i]
				);
			}
		}	
		
		// Draw Background behind labels
		outerRadius = renderBackground();
		
		// Draw Ribbon track highlights
		if(Smits.PhyloCanvas.Render.Parameters.binaryCharts.length){
			var binaryCharts = Smits.PhyloCanvas.Render.Parameters.binaryCharts;
			for(var i in binaryCharts){
				outerRadius = renderBinaryChart(outerRadius, binaryCharts[i].chart, binaryCharts[i]);
			}
		}

		// Draw Bar Chart
		if(Smits.PhyloCanvas.Render.Parameters.barCharts.length){
			var barCharts = Smits.PhyloCanvas.Render.Parameters.barCharts;
			for(var i in barCharts){
				outerRadius = renderBarChart(outerRadius, barCharts[i].chart, barCharts[i]);
			}
		}		

	}
})();

Smits.PhyloCanvas.Render.CircularPhylogram.prototype = {

};/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
var XMLObjectifier = (function() {
	var _clone = function(obj){
		if(!!obj && typeof(obj)==="object"){
			function F(){}
			F.prototype = obj;
			return new F();
		}		
	};
	//Is Numeric check
	var isNumeric = function(s) {
		var testStr = "";
		if(!!s && typeof(s) === "string") { testStr = s; }
		var pattern = /^((-)?([0-9]*)((\.{0,1})([0-9]+))?$)/;
		return pattern.test(testStr);
	};
	var _self = {
	xmlToJSON: function(xdoc) {
		try {
			if(!xdoc){ return null; }
			var tmpObj = {};
				tmpObj.typeOf = "JSXBObject";
			var xroot = (xdoc.nodeType == 9)?xdoc.documentElement:xdoc;
				tmpObj.RootName = xroot.nodeName || "";
			if(xdoc.nodeType == 3 || xdoc.nodeType == 4) {
				return xdoc.nodeValue;
			}
			//Trim function
			function trim(s) {
				return s.replace(/^\s+|\s+$/gm,'');
			}						
			//Alters attribute and collection names to comply with JS
			function formatName(name) {
				var regEx = /-/g;
				var tName = String(name).replace(regEx,"_");
				return tName;
			}
			//Set Attributes of an object
			function setAttributes(obj, node) {
				if(node.attributes.length > 0) {
					var a = node.attributes.length-1;
					var attName;
					obj._attributes = [];
					do { //Order is irrelevant (speed-up)
						attName = String(formatName(node.attributes[a].name));
						obj._attributes.push(attName);				
						obj[attName] = trim(node.attributes[a].value);
					} while(a--);
				}
			}
			
			//Node Prototype
			var _node = (function() {
					var _self = {
						activate: function() {
							var nodes = [];
							if(!!nodes) {
									nodes.getNodesByAttribute = function(attr, obj) {
										if(!!nodes && nodes.length > 0) {
											var out = [];
											var cNode;
											var maxLen = nodes.length -1;
											try {
												do {
													cNode = nodes[maxLen];
													if(cNode[attr] === obj) {
														out.push(cNode);
													}
												} while(maxLen--);
												out.reverse();
												return out;
											} catch(e) {return null;}
											return null;
										}
									};
									nodes.getNodeByAttribute = function(attr, obj) {
										if(!!nodes && nodes.length > 0) {
											var cNode;
											var maxLen = nodes.length -1;
											try {
												do {
													cNode = nodes[maxLen];
													if(cNode[attr] === obj) {
														return cNode;
													}
												} while(maxLen--);
											} catch(e) {return null;}
											return null;
										}
									};
									nodes.getNodesByValue = function(obj) {
										if(!!nodes && nodes.length > 0) {
											var out = [];
											var cNode;
											var maxLen = nodes.length -1;
											try {
												do {
													cNode = nodes[maxLen];
													if(!!cNode.Text && cNode.Text === obj) {
														out.push(cNode);
													}
												} while(maxLen--);
												return out;
											} catch(e) {return null;}
											return null;
										}
									};
									nodes.contains = function(attr, obj) {
										if(!!nodes && nodes.length > 0) {
											var maxLen = nodes.length -1;
											try {
												do {
													if(nodes[maxLen][attr] === obj) {
														return true;
													}
												} while(maxLen--);
											} catch(e) {return false;}
											return false;
										}
									};
									nodes.indexOf = function(attr, obj) {
										var pos = -1;
										if(!!nodes && nodes.length > 0) {
											var maxLen = nodes.length -1;
											try {
												do {
													if(nodes[maxLen][attr] === obj) {
														pos = maxLen;
													}
												} while(maxLen--);
											} catch(e) {return -1;}
											return pos;
										}
									};
									nodes.SortByAttribute = function(col, dir) {
										if(!!nodes && nodes.length > 0) {				
											function getValue(pair, idx) {
												var out = pair[idx];
												out = (bam.validation.isNumeric(out))?parseFloat(out):out;
												return out;
											}
											function sortFn(a, b) {
												var tA, tB;
												tA = getValue(a, col);
												tB = getValue(b, col);
												var res = (tA<tB)?-1:(tB<tA)?1:0;
												if(!!dir) {
													res = (dir.toUpperCase() === "DESC")?(0 - res):res;
												}
												return res;
											}
											nodes.sort(sortFn);
										}
									};
									nodes.SortByValue = function(dir) {
										if(!!nodes && nodes.length > 0) {
											function getValue(pair) {
												var out = pair.Text;
												out = (bam.validation.isNumeric(out))?parseFloat(out):out;
												return out;
											}
											function sortFn(a, b) {
												var tA, tB;
												tA = getValue(a);
												tB = getValue(b);
												var res = (tA<tB)?-1:(tB<tA)?1:0;
												if(!!dir) {
													res = (dir.toUpperCase() === "DESC")?(0 - res):res;
												}
												return res;
											}
											nodes.sort(sortFn);
										}
									};
									nodes.SortByNode = function(node, dir) {
										if(!!nodes && nodes.length > 0) {
											function getValue(pair, node) {
												var out = pair[node][0].Text;
												out = (bam.validation.isNumeric(out))?parseFloat(out):out;
												return out;
											}
											function sortFn(a, b) {										
												var tA, tB;
												tA = getValue(a, node);
												tB = getValue(b, node);
												var res = (tA<tB)?-1:(tB<tA)?1:0;
												if(!!dir) {
													res = (dir.toUpperCase() === "DESC")?(0 - res):res;
												}
												return res;
											}
											nodes.sort(sortFn);
										}
								  };
							}
							return nodes;
						}
					};
					return _self;
			})();
			//Makes a new node of type _node;
			var makeNode = function() {
				var _fn = _clone(_node);					
				return _fn.activate();
			}
			//Set collections
			function setHelpers(grpObj) {
				//Selects a node withing array where attribute = value
				grpObj.getNodeByAttribute = function(attr, obj) {
					if(this.length > 0) {
						var cNode;
						var maxLen = this.length -1;
						try {
							do {
								cNode = this[maxLen];
								if(cNode[attr] == obj) {
									return cNode;
								}
							} while(maxLen--);
						} catch(e) {return false;}
						return false;
					}
				};
				
				grpObj.contains = function(attr, obj) {
					if(this.length > 0) {
						var maxLen = this.length -1;
						try {
							do {
								if(this[maxLen][attr] == obj) {
									return true;
								}
							} while(maxLen--);
						} catch(e) {return false;}
						return false;
					}
				};
				
				grpObj.indexOf = function(attr, obj) {
					var pos = -1;
					if(this.length > 0) {
						var maxLen = this.length -1;
						try {
							do {
								if(this[maxLen][attr] == obj) {
									pos = maxLen;
								}
							} while(maxLen--);
						} catch(e) {return -1;}
						return pos;
					}
				};
				
				grpObj.SortByAttribute = function(col, dir) {
					if(this.length) {				
						function getValue(pair, idx) {
							var out = pair[idx];
							out = (isNumeric(out))?parseFloat(out):out;
							return out;
						}
						function sortFn(a, b) {
							var res = 0;
							var tA, tB;						
							tA = getValue(a, col);
							tB = getValue(b, col);
							if(tA < tB) { res = -1;	} else if(tB < tA) { res = 1; }
							if(dir) {
								res = (dir.toUpperCase() == "DESC")?(0 - res):res;
							}
							return res;
						}
						this.sort(sortFn);
					}
				};
				
				grpObj.SortByValue = function(dir) {
					if(this.length) {
						function getValue(pair) {
							var out = pair.Text;
							out = (isNumeric(out))?parseFloat(out):out;
							return out;
						}
						function sortFn(a, b) {
							var res = 0;
							var tA, tB;
							tA = getValue(a);
							tB = getValue(b);
							if(tA < tB) { res = -1;	} else if(tB < tA) { res = 1; }
							if(dir) {
								res = (dir.toUpperCase() == "DESC")?(0 - res):res;
							}
							return res;
						}
						this.sort(sortFn);
					}
				};
				
				grpObj.SortByNode = function(node, dir) {
					if(this.length) {
						function getValue(pair, node) {
							var out = pair[node][0].Text;
							out = (isNumeric(out))?parseFloat(out):out;
							return out;
						}
						function sortFn(a, b) {
							var res = 0;
							var tA, tB;
							tA = getValue(a, node);
							tB = getValue(b, node);
							if(tA < tB) { res = -1;	} else if(tB < tA) { res = 1; }
							if(dir) {
								res = (dir.toUpperCase() == "DESC")?(0 - res):res;
							}
							return res;
						}
						this.sort(sortFn);
					}
				};
			}
			//Recursive JSON Assembler
			//Set Object Nodes
			function setObjects(obj, node) {
				var elemName;	//Element name
				var cnode;	//Current Node
				var tObj;	//New subnode
				var cName = "";
				if(!node) { return null; }				
				//Set node attributes if any
				if(node.attributes.length > 0){setAttributes(obj, node);}				
				obj.Text = "";
				if(node.hasChildNodes()) {
					var nodeCount = node.childNodes.length - 1;	
					var n = 0;
					do { //Order is irrelevant (speed-up)
						cnode = node.childNodes[n];
						switch(cnode.nodeType) {
							case 1: //Node
							//Process child nodes
							obj._children = [];
							//SOAP XML FIX to remove namespaces (i.e. soapenv:)
							elemName = (cnode.localName)?cnode.localName:cnode.baseName;
							elemName = formatName(elemName);
							if(cName != elemName) { obj._children.push(elemName); }
								//Create sub elemns array
								if(!obj[elemName]) {
									obj[elemName] = []; //Create Collection
								}
								tObj = {};
								obj[elemName].push(tObj);
								if(cnode.attributes.length > 0) {
									setAttributes(tObj, cnode);
								}
								//Set Helper functions (contains, indexOf, sort, etc);
								if(!obj[elemName].contains) {
									setHelpers(obj[elemName]);
								}	
							cName = elemName;
							if(cnode.hasChildNodes()) {
								setObjects(tObj, cnode); //Recursive Call
							}
							break;
							case 3: //Text Value
							obj.Text += trim(cnode.nodeValue);
							break;
							case 4: //CDATA
							obj.Text += (cnode.text)?trim(cnode.text):trim(cnode.nodeValue);
							break;
						}
					} while(n++ < nodeCount);
				}
			}			
			//RUN
			setObjects(tmpObj, xroot);
			//Clean-up memmory
			xdoc = null;
			xroot = null;
			return tmpObj;	
		} catch(e) {
				return null;	
		}	
	},

	//Converts Text to XML DOM
	textToXML: function(strXML) {
		var xmlDoc = null;
		try {
			xmlDoc = (document.all)?new ActiveXObject("Microsoft.XMLDOM"):new DOMParser();
			xmlDoc.async = false;
		} catch(e) {throw new Error("XML Parser could not be instantiated");}
		var out;
		try {
			if(document.all) {
				out = (xmlDoc.loadXML(strXML))?xmlDoc:false;
			} else {		
				out = xmlDoc.parseFromString(strXML, "text/xml");
			}
		} catch(e) { throw new Error("Error parsing XML string"); }
		return out;
	}
	};
	return _self;
})();
