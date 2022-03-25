from neo4j import GraphDatabase
from json import dumps
import re
import logging

class GraphAPI:

    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth = (user, password))

    def close(self):
        self.driver.close()

    def writeGraph(self, relationships, node_properties, additional_nodes):
        with self.driver.session() as session:
            for node in additional_nodes:
                result = session.write_transaction(self.AddNode, 
                                                   node, 
                                                   node_properties)
            for triplet in relationships:
                result = session.write_transaction(self.AddRelationship, 
                                                   triplet, 
                                                   node_properties)

    @staticmethod
    def DictToString(d):
        # convert dictionary to string
        d = dumps(d)
        # remove quotes around dictionary keys
        d = re.sub('"(?=:)|(?<={)"|(?<=, )"', '', d)
        # add whitespace after brackets
        d = re.sub('({)(\w)', r'\1 \2', d)
        d = re.sub('(\")(})', r'\1 \2', d)
        return d
              
    @staticmethod
    ## writes a cypher query to add a relationship to the graph - assumes link comes as a dictionary
    def AddRelationship(tx, triplet, node_properties):
        
        # split triplet into head, link, tail
        head, link, tail = triplet
        
        # get head and tail properties from the properties dictionary, convert to str 
        head_properties = GraphAPI.DictToString(node_properties[head]['properties'])
        tail_properties = GraphAPI.DictToString(node_properties[tail]['properties'])
        
        # get node type for head and tail
        head_type = node_properties[head]['type']
        tail_type = node_properties[tail]['type']
        
        # write query using head, link, and tail properties 
        link_type = link['type']
        link_properties = GraphAPI.DictToString(link['properties'])
        query = (
            f"MERGE (h: {head_type} {head_properties}) " 
            f"MERGE (t: {tail_type} {tail_properties}) "
            f"MERGE (h)-[l:{link_type} {link_properties}]->(t) "
            "return h, l, t"
        )
        result = tx.run(query)    
        
        try:
            return [{'head': record['h'], 'tail': record['t'], 'link': record['l']}
                    for record in result]
        # Capture any errors along with the query and data for traceability
        except ServiceUnavailable as exception:
            logging.error("{query} raised an error: \n {exception}".format(
                query=query, exception=exception))
            raise
            
    @staticmethod
    def AddNode(tx, node, node_properties):
        
        props = GraphAPI.DictToString(node_properties[node]['properties'])
        node_type = node_properties[node]['type']
        
        query = (
            f"MERGE (e: {node_type} {props}) " 
            "return e"
        )
        result = tx.run(query)
        
        try:
            return [{'node': record['e']} for record in result]
        # Capture any errors along with the query and data for traceability
        except ServiceUnavailable as exception:
            logging.error("{query} raised an error: \n {exception}".format(
                query=query, exception=exception))
            raise