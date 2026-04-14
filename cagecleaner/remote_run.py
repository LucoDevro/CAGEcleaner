#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cagecleaner.run import Run

import logging
from abc import abstractmethod


LOG = logging.getLogger()


class RemoteRun(Run):
    """
    Abstract intermediary class grouping the methods shared by every run involving remote sequence files.
    
    Inherits from:
        Run: Base class providing argument parsing, hit recovery, session filtering and output generation functionalities
        
     See Also:
         RemoteRegionRun: Region-based dereplication for hits in remote sequences.
         RemoteGenomeRun: Whole-genome dereplication for hits in remote sequences.
    """
    
    def __init__(self, args):
        """
        Initialise a RemoteRun instance.
        
        Runs the base class init.
        Excludes scaffolds from the analysis as specified by the user.
        
        Args:
            args (argparse.Namespace): Parsed command-line arguments
            
        Returns:
            None
        """
        # Call the parent constructor:
        super().__init__(args)
        
        # Remove Organisms specified by the user:
        if self.excluded_organisms != {''}:
            # Replace colons with spaces to get the correct matching in the binary table.
            self.excluded_organisms = {org.replace(':', ' ') for org in self.excluded_organisms}
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
    def join_dereplication_with_binary(self):
        """
        Join dereplication results with the binary table.
        
        Mutates:
            self.binary_df: Adds 'representative' and 'dereplication_status' columns.
            
         Expected Result:
             self.binary_df should now have columns:
             - representative: Genome ID of the dereplication representative
             - dereplication_status: 'dereplication_representative' | 'redundant'
             
         Returns:
             None
             
        Notes:
            This is the abstract method inherited from the Run parent class.
            It is not meant to be implemented at this level. Only the child classes inheriting this method
            are expected to provide a workflow-dependent specific implementation.
        """
        pass
    
    