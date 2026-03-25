#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run

import logging
from abc import abstractmethod


LOG = logging.getLogger()


class RemoteRun(Run):
    
    def __init__(self, args):
        # Call the parent constructor:
        super().__init__(args)
        
        # Remove Organisms specified by the user:
        if self.excluded_organisms != {''}:
            # Replace colons with spaces to get the correct matching in the binary table.
            self.excluded_organisms = {org.replace(':', ' ', regex=False) for org in self.excluded_organisms}
            LOG.debug(f"Excluding the following organisms: {', '.join(self.excluded_organisms)}")
            # Exclude them:
            self.binary_df = self.binary_df[~self.binary_df['Organism'].isin(self.excluded_organisms)]
        
        # Remove scaffolds that the user wants excluded:
        if self.excluded_scaffolds != {''}:
            LOG.debug(f"Excluding the following scaffolds: {', '.join(self.excluded_scaffolds)}")
            self.binary_df = self.binary_df[~self.binary_df['Scaffold'].isin(self.excluded_scaffolds)]
        
        # Replace colons in the bypass assemblies as wel:
        if self.bypass_organisms != {''}:
            self.bypass_organisms = {org.replace(':', ' ') for org in self.bypass_organisms}
            
        return None
    
    @abstractmethod
    def mapDereplicationToBinary(self):
        pass
    
    