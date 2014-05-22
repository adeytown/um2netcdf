/**============================================================================
                 U M 2 N e t C D F  V e r s i o n 2 . 0
                 --------------------------------------

    Main author: Mark Cheeseman
                 National Institute of Water & Atmospheric Research (Ltd)
                 Wellington, New Zealand
                 February 2014

    UM2NetCDF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    UM2NetCDF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the main UM2NetCDF
    directory.  Alternatively, please see <http://www.gnu.org/licenses/>.
 **============================================================================*/

#include <string.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "field_def.h"

/***
 *** PARSE_ITEM 
 ***
 *** Subroutine that reads all the metadata characteristics of a variable 
 *** defined in the user-specified stash file.  
 ***
 ***   Mark Cheeseman, NIWA
 ***   November 29, 2013
 ***/

void parse_item( xmlDocPtr doc, xmlNodePtr cur, um_field_metadata *fd, int cnt ) {

     xmlChar *str;
     cur = cur->xmlChildrenNode;
     while ( cur!=NULL ) {
          if ((!xmlStrcmp(cur->name, (const xmlChar *)"stash_code"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             fd[cnt].code = atoi((const char *) str);
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"varname"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) {
                memset( fd[cnt].varname, '\0', sizeof(fd[cnt].varname) );
                strcpy( fd[cnt].varname,(const char *)str );
             }
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"longname"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) { 
                memset( fd[cnt].longname, '\0', sizeof(fd[cnt].longname) );
                strcpy( fd[cnt].longname,(const char *)str ); 
             }
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"standardname"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) { 
                memset( fd[cnt].stdname, '\0', sizeof(fd[cnt].stdname) );
                strcpy( fd[cnt].stdname,(const char *)str ); 
             }
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"units"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) { 
                memset( fd[cnt].units, '\0', sizeof(fd[cnt].units) );
                strcpy( fd[cnt].units,(const char *)str ); 
             }
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"validmax"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) { fd[cnt].validmax = atof((const char *)str); }
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"validmin"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) { fd[cnt].validmin = atof((const char *)str); }
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"scalefact"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) { fd[cnt].scale = atof((const char *)str); }
             else             { fd[cnt].scale = 1.0; }
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"umgrid"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             fd[cnt].umgrid = atoi((const char *) str);
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"accum_field"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             fd[cnt].accum = atoi((const char *) str);
             xmlFree( str );
          } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"level_type"))) {
             str = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
             if ( str!=NULL ) { fd[cnt].level_type = atoi((const char *) str); } 
             else             { fd[cnt].level_type = 2; } 
             xmlFree( str );
          }
          cur = cur->next;
     }
     return;

}

/***
 *** READ_STASH_FILE 
 ***
 *** Subroutine that opens the user-specified XML stash file and reads 
 *** in the metadata associated with each variable defined in the file. 
 ***
 ***   Mark Cheeseman, NIWA
 ***   November 29, 2013
 ***/

int read_stash_file( char *filename ) {

     int cnt, model_num, section_num;
     xmlDocPtr doc;
     xmlNodePtr model, section, item;
     xmlChar *id;

/**
 ** Open the XML file and verify its layout
 **---------------------------------------------------------------------------*/
     doc = xmlParseFile( filename );
     if ( doc==NULL ) { return 0; }
     else {
          model = xmlDocGetRootElement(doc);
          if (model == NULL) { return 0; }
          else { 
              if (xmlStrcmp(model->name, (const xmlChar *) "stash2cf")) { return 0; } 
          }
     }           

/**
 ** Determine number of variables defined in XML file (num_xml_vars) 
 **---------------------------------------------------------------------------*/
     model = model->xmlChildrenNode;
     num_xml_vars = 0;
     while ( model!= NULL ) {
           if ((!xmlStrcmp(model->name, (const xmlChar *)"model"))) {
              section = model->xmlChildrenNode;
              while ( section!= NULL ) {
                    if ((!xmlStrcmp(section->name, (const xmlChar *)"section"))) {
                       item = section->xmlChildrenNode;
                       while ( item!= NULL ) {
                             if ((!xmlStrcmp(item->name, (const xmlChar *)"item"))) { num_xml_vars++; }
                             item = item->next;
                       }
                    }
                     section = section->next;
              }
           }
           model = model->next;
     }
     xmlFreeDoc( doc );

/**
 ** Reopen XML file 
 **---------------------------------------------------------------------------*/
     doc = xmlParseFile( filename );
     model = xmlDocGetRootElement(doc);

/**
 ** Allocate memory to hold metadata for each defined variable  
 **---------------------------------------------------------------------------*/
     um_vars = (um_field_metadata *)malloc( num_xml_vars*sizeof(um_field_metadata) );

/**
 ** Loop over the different models
 **---------------------------------------------------------------------------*/
     cnt = 0;
     model = xmlDocGetRootElement(doc);
     model = model->xmlChildrenNode;
     while ( model!= NULL ) {
           if ((!xmlStrcmp(model->name, (const xmlChar *)"model"))) {
              id = xmlGetProp( model, (const xmlChar *)"model_id" );
              model_num = atoi((const char *)id);
              xmlFree( id );

/**
 ** Now loop over the sections within that model's domain
 **---------------------------------------------------------------------------*/
              section = model->xmlChildrenNode;
              while ( section!= NULL ) {
                    if ((!xmlStrcmp(section->name, (const xmlChar *)"section"))) {
                       id = xmlGetProp( section, (const xmlChar *)"section_id" );
                       section_num = atoi((const char *)id);
                       xmlFree( id );

/**
 ** Now loop over the individual items within each section
 **---------------------------------------------------------------------------*/
                       item = section->xmlChildrenNode;
                       while ( item!= NULL ) {
                             if ((!xmlStrcmp(item->name, (const xmlChar *)"item"))) {
                                um_vars[cnt].model = model_num;
                                um_vars[cnt].section = section_num;
                                parse_item( doc, item, um_vars, cnt );
                                cnt++;
                             }
                             item = item->next;
                       } 

                    }
                    section = section->next;
              } 
           }
           model = model->next;
     }
     xmlFreeDoc( doc );

     return 1; 
}


/***
 *** READ_CONFIG_FILE
 ***
 *** Subroutine that opens the user-specified XML config file and reads
 *** in the metadata associated with the UM run.
 ***
 ***   Mark Cheeseman, NIWA
 ***   May 22, 2014
 ***/

int read_config_file( char *filename ) {

     xmlDocPtr doc;
     xmlNodePtr run;
     xmlChar *str;

/**
 ** Open the XML file
 **---------------------------------------------------------------------------*/
     doc = xmlParseFile( filename );
     if ( doc==NULL ) { return 0; }

     run = xmlDocGetRootElement(doc);
     if (run == NULL) { return 0; }
     else { if (xmlStrcmp(run->name, (const xmlChar *) "run_config")) { return 0; } }

/**
 ** Loop over the different models
 **---------------------------------------------------------------------------*/
     run = run->xmlChildrenNode;
     while ( run!= NULL ) {
           if ((!xmlStrcmp(run->name, (const xmlChar *)"institution"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) {
                 memset( run_config.institution, '\0', sizeof(run_config.institution) );
                 strcpy( run_config.institution,(const char *)str );
              }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"ps"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) { run_config.eps = atoi((const char *) str); }
              else             { run_config.eps = 34; }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"niwa_eps"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) { run_config.eps = atoi((const char *) str); }
              else             { run_config.eps = 0; }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"rose_id"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) { run_config.rose_id = atoi((const char *) str); }
              else             { run_config.rose_id = 0; }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"model_name"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) {
                 memset( run_config.model, '\0', sizeof(run_config.model) );
                 strcpy( run_config.model,(const char *)str );
              }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"references"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) {
                 memset( run_config.ref, '\0', sizeof(run_config.ref) );
                 strcpy( run_config.ref,(const char *)str );
              }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"comment"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) {
                 memset( run_config.comment, '\0', sizeof(run_config.comment) );
                 strcpy( run_config.comment,(const char *)str );
              }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"title"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) {
                 memset( run_config.title, '\0', sizeof(run_config.title) );
                 strcpy( run_config.title,(const char *)str );
              }
              xmlFree( str );
           }
           if ((!xmlStrcmp(run->name, (const xmlChar *)"data_assimilation_method"))) {
              str = xmlNodeListGetString(doc, run->xmlChildrenNode, 1);
              if ( str!=NULL ) {
                 memset( run_config.assim, '\0', sizeof(run_config.assim) );
                 strcpy( run_config.assim,(const char *)str );
              }
              xmlFree( str );
           }
           run = run->next;
     }
     xmlFreeDoc( doc );

     return 1;
}


