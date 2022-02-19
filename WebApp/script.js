function post_cypherquery(cypher_query,update){
  // while busy, show we're doing something in the messageArea.
  //  $('#messageArea').html('<h3>(loading)</h3>');


   $.ajaxSetup({
  headers: {
    "Authorization": 'Basic ' + window.btoa("neo4j"+":"+"Cellular senescence")
  }
  })
  // get the data from neo4j
   $.ajax({
       url: "http://ss09.hutchison-mrc.cam.ac.uk:7474/db/data/transaction/commit",
       type: 'POST',
       data: JSON.stringify({ "statements": [{ "statement": cypher_query }] }),
       contentType: 'application/json',
       accept: 'application/json; charset=UTF-8',
       success: function () { },
       error: function (jqXHR, textStatus, errorThrown) { $('#messageArea').html('<h3>' + textStatus + ' : ' + errorThrown + '</h3>') },
       complete: function () { }
   })
   .then(function (data) {
    //  $('#outputArea').html("<p>Query: '"+ cypher_query +"'</p>");
    //  $('#messageArea').html('');

// save the results from neo4j in the d3_data object
     var d3_data = [];
     $.each(data.results[0].data, function(k,v){d3_data.push(v.row);});

// link_data has the data structure which can be read by d3.js
    var link_data = d3_data[0][0];
    console.log(link_data);

// If the result is empty, alert the message
    if(link_data.nodes.length == 0){
      alert("Sorry, no results were found for your query!");
    }

// where the second parameter is passed to the function, update the attribute of nodes for visualizing them in different colors later
    if(update == "update_yes"){
      update_color(link_data);
      // console.log(link_data);
    };


// when the window is changed, reload the graph accordingly
d3.select(window).on("resize",callFunction);
callFunction();

function callFunction(){

  var svgtest = d3.select("body").select("svg");
  if (!svgtest.empty()){
    svgtest.remove();
  }

// define the size of the svg object on which the nodes and lines will be attached on later
    var width = window.innerWidth - document.getElementById('sidebar').clientWidth;
    var height = window.innerHeight-100;
    var radius = 25;
    var margin = {top:40, right:35, bottom:20, left:10};

    var graph = d3.select('#outputArea')
                .append("svg")
                  .attr("class","background")
                  .attr("width",width - margin.left - margin.right)
                  .attr("height", height)
                  ;
                // .append("g")
                //   .attr("transform", "translate(" + margin.left + "," + margin.right +")");

    // var color = d3.scaleOrdinal(d3.schemeCategory10);
    // console.log(color(1));

    var color = d3.scaleOrdinal(d3.schemeCategory20b);
    console.log(color(1));

    var color_array = ["#3d090a","#ba282b","#306649","#3a258c","#7c420e"];
    function color2(i){
      return color_array[i];
    }
    console.log(color2(1));


    var simulation = d3.forceSimulation()
    .force("link", d3.forceLink().id(function(d) { return d.id; }))
    .force("charge", d3.forceManyBody())
    // .strength(-10))
    .force("center", d3.forceCenter(width / 2, height / 2))
    .force("collision", d3.forceCollide().radius(13))

    var zoom_handler = d3.zoom()
                .scaleExtent([1,3])
                .on("zoom",zoomed);

    var drag_handler = d3.drag()
      .on("start", dragstarted)
      .on("drag", dragged)
      .on("end", dragended);

// Define the variables of links, nodes, and labels of the nodes to be attached to the svg
    var link = graph.selectAll(".links");
    var node = graph.selectAll(".nodes");
    var nodelabels = graph.selectAll(".nodelabel");

    graph.call(zoom_handler).on("dblclick.zoom",null);

    graph.append("title")
        .text("Genes");

    graph.append("desc")
        .text("This is a demo to buidling visualization using d3 for neo4j data");

// Bind the array "links" to d3 lines
    link = link.data(link_data.links).enter()
                .append("line")
                .attr("class","links")
                .attr("stroke",function(d){return color(d.type);})
                .attr("stroke-opacity","0.7")
                .attr("stroke-width","1.5")
            .on("mouseover", function(d){
                d3.select(this)
                    .style("stroke-width","3")
              })
            .on("mouseout", function(d){
              d3.selectAll("line")
              .style("stroke-width","1.5")
            });

// Bind the array "nodes" to d3 circles
    node = node.data(link_data.nodes)
            .enter()
            .append("circle")
            .attr("class","nodes")
            // .style("fill","#6d7fcc")
            .style("fill",function(d){
                return "#d62728";
            })
            .style("opacity",0.7)
            .attr("stroke",function(d){
                return "#037013";
            })
            .attr("stroke-width",1.5)
            .attr("r",10)
            .on("dblclick",releasenode)
            .call(drag_handler);

// Bind the array "nodes" to d3 text as the labels for the nodes
    nodelabels = nodelabels.data(link_data.nodes)
            .enter()
            .append("text")
                 .attr("class","nodelabel")
                 .attr("visibility","hidden")
                 .attr("font-family","sans-serif")
                 .attr("font-size","10px")
                //  .text(function(d.id);)
                 .text(function(d){
                   return d.symbol;
                 });

// Append the senescence types and logFC values to the nodes
     node.append("title")
         .text(function(d){
           if(d.differential_expression != "None"){
             return "Differentially-expressed in " + d.differential_expression + " senescence";
           }else{
             return "Not differentially-expressed in senescence, sourced from databases " + d.source;
           }
         })

// Append the relation type to the lines as titles
     link.append("title")
         .text(function(d){
           return d.type;
         })

// Show and hide the lables for the nodes when clicking anywhere on the graph
    graph.on('click', change_label_visibility);
// When double clicking anywhere on the svg, the graph will be resetted to its original size and position
    graph.on('dblclick',zoom_reset);
    function zoom_reset(){
      graph.transition()
       .duration(750)
       .call(zoom_handler.transform, d3.zoomIdentity);
    }


      // console.log(nodelabels.attr("visibility"));


     simulation.nodes(link_data.nodes)
              .on("tick", ticked);

     simulation.force("link")
              .links(link_data.links);


      function ticked() {

        node
        // contrain the area where the nodes are supposed to be shown in
            .attr("cx", function(d) { return d.x = Math.max(radius + 2*margin.left, Math.min(width - radius - 4*margin.right, d.x));  })
            .attr("cy", function(d) { return d.y = Math.max(radius + margin.top, Math.min(height - radius - margin.bottom, d.y)); });

        link
            .attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
            .attr("x2", function(d) { return d.target.x; })
            .attr("y2", function(d) { return d.target.y; });

        nodelabels
                  .attr("x", function(d) { return d.x + 8; })
                  .attr("y", function(d) { return d.y; })
                  .attr("class","noselect");
                }


          function dragstarted(d) {
              if (!d3.event.active) simulation.alphaTarget(0.3).restart();
              d.fx = d.x;
              d.fy = d.y;
            }

          function dragged(d) {
            d.fx = d3.event.x;
            d.fy = d3.event.y;
          }

          function dragended(d) {
            if (!d3.event.active) simulation.alphaTarget(0.0);
            d.fx = d3.event.x;
            d.fy = d3.event.y;
          }

          function releasenode(d){
            d.fx = null;
            d.fy = null;
          }

          function zoomed(){
            graph.attr("transform", d3.event.transform);
          }

// function for changing the visibility of the node labels
          function change_label_visibility(){
            if (nodelabels.attr("visibility") == "visible"){
              nodelabels.attr("visibility","hidden");
            }
            else{
              nodelabels.attr("visibility","visible");
            }
          }

        }//end of resize
      });//end of data-node binding

  };//end of post_cypherquery


// Function for get the values of the selected options and save them in an array
function getCheckedValues_1_1(){
  var chkArray = [];


  $("#checkboxlist_1_1 input:checked").each(function(){
    chkArray.push(":" + $(this).val());
  });

  var selected_types;
  selected_types = chkArray.join('|');
  console.log(selected_types);
  return(selected_types);
};


getRadioButtonValue1();
function getRadioButtonValue1(){
  var type_selected = $('input[name="result1"]:checked').val();
  console.log(type_selected);
  return(type_selected);
};

// Function for searching genes in the Neo4j database using a list of gene symbols
function search_gene(){
  var inputstring = document.getElementById('nodename').value;
  var inputarray = inputstring.split(",");
  var input = "[";
  for(i =0 ; i < inputarray.length-1; i++){
    input += '"' + inputarray[i] + '",';
    // inputarray[i] = '"' + inputarray[i] + '"' ;
  };
  input += '"' + inputarray[inputarray.length-1] + '"';
  input += "]";

  var selected = "[" + getCheckedValues_1_1() + "]";
  // console.log(selected);

  var type_selected = getRadioButtonValue1();

  if(selected.length <= 2){
    alert("Please choose relation types");
  }else {
    var search_command = "MATCH path = (a:Gene)-"+selected+"->(b:Gene) WHERE a.symbol in "+input+ " UNWIND nodes(path) as p UNWIND relationships(path) as r RETURN {nodes: collect(DISTINCT p{.*, id: id(p), label:labels(p), oncogene_LFC: p.oncogene_LFC, DNA_damage_LFC: p.DNAdamage_LFC, replicative_LFC: p.replicative_LFC}), links: collect(DISTINCT {source: id(startNode(r)), target: id(endNode(r)), type: type(r)})}"
    post_cypherquery(search_command);
  }
};


// Function for parsing the input text and returning an array containing gene symbols
function get_gene_list(){
  var raw_input = document.getElementById('nodename_plus').value;
  // console.log(raw_input);
  var inputarray = raw_input.split("\n");
  // console.log(inputarray);
  var input_gene=[];
  for (var i = 0; i < inputarray.length; i++){
    var input_record = inputarray[i];
    input_record_splitted = input_record.split(",");
    input_gene.push(input_record_splitted[0]);
  }
  // console.log(input_name);
  return input_gene;
}

// Function for parsing the input text and returning an array containing log fold change values
function get_logFC_list(){
  var raw_input = document.getElementById('nodename_plus').value;
  var inputarray = raw_input.split("\n");
  var input_logFC=[];
  for (var i = 0; i < inputarray.length; i++){
    var input_record = inputarray[i];
    input_record_splitted = input_record.split(",");
    input_logFC.push(input_record_splitted[1]);
  }
  return input_logFC;
}


// Function for get the values of the selected options in select box (2_1) and save them in an array
function getCheckedValues_2_1(){
  var chkArray = [];

  $("#checkboxlist_2_1 input:checked").each(function(){
    chkArray.push(":" + $(this).val());
  });

  var selected_types;
  selected_types = chkArray.join('|');
  console.log(selected_types);
  return(selected_types);
};


function getRadioButtonValue2(){
  var type_selected = $('input[name="result"]:checked').val();
  console.log(type_selected);
  return(type_selected);
};


//funciton for searching genes (highlighted according to logFC)
function search_with_highlight(){

  var input_gene = get_gene_list();
  var input = "[";
  for(i =0 ; i < input_gene.length-1; i++){
    input += '"' + input_gene[i] + '",';
    // inputarray[i] = '"' + inputarray[i] + '"' ;
  };
  input += '"' + input_gene[input_gene.length-1] + '"';
  input += "]";
  console.log(input);

  var selected = "[" + getCheckedValues_2_1() + "]";
  console.log(selected);

  var type_selected = getRadioButtonValue2();


  if(selected.length <= 2){
    alert("Please choose relation types");
  }else {
    if(type_selected == "non-target"){
      var search_command = "MATCH path = (a:Gene)-"+selected+"->(b:Gene) WHERE a.symbol in "+input+ " AND size(b.IFN_type)=0 UNWIND nodes(path) as p UNWIND relationships(path) as r RETURN {nodes: collect(DISTINCT p{.*, id: id(p), label:labels(p), logFC:null}), links: collect(DISTINCT {source: id(startNode(r)), target: id(endNode(r)), type: type(r)})}"
    }else if(type_selected == "targets"){
      var search_command = "MATCH path = (a:Gene)-"+selected+"->(b:Gene) WHERE a.symbol in "+input+ " AND size(b.IFN_type)<>0 UNWIND nodes(path) as p UNWIND relationships(path) as r RETURN {nodes: collect(DISTINCT p{.*, id: id(p), label:labels(p), logFC:null}), links: collect(DISTINCT {source: id(startNode(r)), target: id(endNode(r)), type: type(r)})}"
    }else{
      var search_command = "MATCH path = (a:Gene)-"+selected+"->(b:Gene) WHERE a.symbol in "+input+ " UNWIND nodes(path) as p UNWIND relationships(path) as r RETURN {nodes: collect(DISTINCT p{.*, id: id(p), label:labels(p), logFC:null}), links: collect(DISTINCT {source: id(startNode(r)), target: id(endNode(r)), type: type(r)})}"
    }
    post_cypherquery(search_command,"update_yes");
  }

};

// Function for updating the logFC attribute the nodes according to the input
function update_color(object){

  var name_list = get_gene_list();
  var color_list = get_logFC_list();
  // console.log(name_list);
  // console.log(object.nodes[1].label[0] == "Movie");

  for (var i = 0; i < name_list.length; i++){
    for (var j = 0; j < object.nodes.length; j++){
        if (name_list[i] == object.nodes[j].symbol){
          object.nodes[j].logFC = color_list[i];
        }
    ;}
  ;}
};

//function for changing visibility of elements
function change_variable_visibility(id){
  var x = document.getElementById(id);
  if (x.style.visibility === 'visible'){
    x.style.visibility = 'hidden';
  }
  else{
    x.style.visibility = 'visible';
  }
};

//function for resetting input
function reset_input(id){
  document.getElementById(id).reset();
}

//toggle button
$(document).ready(function(){

  $('#sidebarCollapse').on('click', function(){
    $('#sidebar').toggleClass('active');
  });
});
