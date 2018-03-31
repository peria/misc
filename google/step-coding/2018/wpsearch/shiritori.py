#!/usr/bin/python3

import json


class Shiritori(object):

    def __init__(self, index):
        self.index = index
        self.category = None
        self.initial = None

    def initialize(self, query):
        if self.index.search_category(query):
            self.category = query
            return self.reply('{0} のカテゴリでしりとりしましょう。何から初めますか。'.format(query))
        else:
            return self.reply('{0} のカテゴリは知りません。'.format(query))

    def reply(self, query):
        return json.dumps({
            'textToSpeech': query
        }, indent=2, separators=(',', ': '), ensure_ascii=False)

    @property
    def is_initialized(self):
        return bool(self.category)
