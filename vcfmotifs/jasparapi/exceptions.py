class MatrixFormatError(Exception):
    """Raised when given format is not supported by jaspar API"""

    def __init__(self, given_format, message="is not valid format"):
        self.format = given_format
        self.message = message
        super(Exception, self).__init__(self.message)

    def __str__(self):
        message_out = "{} -> {} \n\t\t supported formats are: jaspar, meme, transfac, pfm".format(self.format,
                                                                                                  self.message)
        return message_out


class MatrixIDError(Exception):
    def __init__(self, matrix_id, message="ID not found"):
        self.matrix_id = matrix_id
        self.message = message
        super(Exception, self).__init__(self.message)

    def __str__(self):
        return "{} : {}".format(self.message, self.matrix_id)


class JasparCollectionError(Exception):
    def __init__(self, collname, message="Collection not found"):
        self.collname = collname
        self.message = message
        super(Exception, self).__init__(self.message)

    def __str__(self):
        valid_colls = "CNE, CORE, FAM, PBM, PBM_HLH, PBM_HOMEO, PHYLOFACTS, POLII, SPLICE, UNVALIDATED"
        return "{} : {} /n/t/t Valid collections are: {}".format(self.message, self.collname, valid_colls)
