function neo4jconnection(search_command, update){

    // login details
    const url = "bolt://ss09.hutchison-mrc.cam.ac.uk:7687";
    const user = "neo4j";
    const password = "doxorubicin";

    const driver = neo4j.driver(url, neo4j.auth.basic(user, password),
                                {
                                    encrypted: 'ENCRYPTION_ON',
                                    disableLosslessIntegers: true,

                                })

    const session = driver.session({

    })

    // example of search term return all the nodes:
    // var search_command = "MATCH (node) return node"

    const result =  session.run(search_command).subscribe({
        onKeys: keys => {
          console.log(keys)
        },
        onNext: record => {
          window.neo4jdata = record["_fields"]
          query(record["_fields"],update)
        },
        onCompleted: () => {
          session.close() // returns a Promise
        },
        onError: error => {
          console.log(error)
        }
      })


    driver.close()

}

function pathwaySelect(pathway){
	console.log('Selecting pathway')
	pathway_selected = $("input[type='radio'][name='pathway_select']:checked").val()
    var search_command = "MATCH path = (a:Gene)-[rel]->(b:Gene) WHERE rel.pathway_name = " + pathway + " return path limit 50"
    neo4jconnection(search_command,"update_yes");
}