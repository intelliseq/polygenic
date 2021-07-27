import json
from typing import List
from typing import Dict


class WebsiteTextData(object):
    def __init__(self, title: str, accessed: str, link: str):
        super().__init__()
        self.title = title.strip('. ')
        self.accessed = accessed.strip('. ')
        self.links = link.strip('. ')


class ArticleLink(object):
    def __init__(self, link: str):
        super().__init__()
        self.link = link.strip('. ')


class References(object):
    def __init__(self, website_text_data_container: List[WebsiteTextData], article_links: List[ArticleLink]):
        super().__init__()
        self.website_text_data = website_text_data_container
        self.article_links = article_links

    def __str__(self):
        return json.dumps(self.to_dict())

    def to_dict(self) -> Dict[str, List]:
        return {k: [sub_val.__dict__ for sub_val in v] for k, v in self.__dict__.items()}


if __name__ == '__main__':
    print(WebsiteTextData('bla', 'acc', 'l1'))
